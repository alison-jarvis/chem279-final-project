#include <iostream> 
#include <string> 
#include <fstream>
#include <sstream>
#include <utility>
#include <armadillo>
#include <cmath>
#include <map>
#include <unordered_map>
#include <vector>
#include <filesystem>
#include <nlohmann/json.hpp> 

namespace fs = std::filesystem;
using json = nlohmann::json; 


/*
DATA STRUCTURES AND LOOKUPS
*/

// Structure for atoms
struct Atom{
    int element;
    arma::vec location;
    double beta;
    int Z;
};

// Structure for primitive gaussian
struct Primitive_Gauss{
    double exp;
    double cc;
};

// Function to read basis from basis path
std::vector<Primitive_Gauss> read_sto3g_basis(std::string sto3g_path){
    // Full basis path, assumes in basis folder, as key_STO3G.json
    std::string full_basis_path = "basis/" + sto3g_path + "_STO3G.json";
    // Read in json file
    std::ifstream basis_file(full_basis_path);
    json basis = json::parse(basis_file);
    // Define empty vector of primitives
    std::vector<Primitive_Gauss> basis_vec;
    // Read the contracted gaussian information out
    for (const auto& prim : basis["contracted_gaussians"]) {
        // Turn this into a primitive gaussian object
        Primitive_Gauss primg;
        primg.exp = prim["exponent"];
        primg.cc = prim["contraction_coefficient"];
        // Add object to basis vector
        basis_vec.push_back(primg);
    }
    return basis_vec;
}

// Struct for atomic sto-3g shell
struct Sto3G_Shell{
    int L;
    double ionization_energy;
    double IA_coeff;
    std::vector<Primitive_Gauss> sto3g;
    Sto3G_Shell(int L, double ie, double iac, std::string sto3g_path) : L(L), ionization_energy(ie), IA_coeff(iac), sto3g(read_sto3g_basis(sto3g_path)){}

};

// Define an object representing an sto-3g basis function
struct Sto3G_Basis{
    int element; // atomic number
    int atomic_index; // atomic index (in atoms vector, for mapping)
    double ionization_energy; // ionization energy
    arma::vec center; // location of atom
    std::vector<int> triplet; // angular momentum triplet
    Sto3G_Shell sto3g; // sto-3g lookup table for basis
};

// Lookup table of atomic number, sto3g shells
std::map<int, std::vector<Sto3G_Shell>> shell_lookup = {
    {1, { {0, -13.6, 7.176, "H_s"} }}, // Hydrogen
    {6, { {0, -21.4, 14.051, "C_s"}, {1, -11.4, 5.572, "C_p"} }}, // Carbon
    {7, { {0, -0.00, 19.316, "N_s"}, {1, -0.00, 7.275, "N_p"} }}, // Nitrogen
    {8, { {0, -0.00, 25.390, "O_s"}, {1, -0.00, 9.111, "O_p"} }}, // Oxygen
    {9, { {0, -0.00, 32.272, "F_s"}, {1, -0.00, 11.080, "F_p"} }} // Fluorine
};

// Lookup table of atomic number, {beta, Z}
std::map<int, std::vector<int>> atom_lookup = {
    {1, {-9, 1}},
    {6, {-21, 4}},
    {7, {-25, 5}},
    {8, {-31, 6}},
    {9, {-39, 7}}
};


/*
HARTREE FOCK CODE
*/

// Parse each xyz file line into element and location
Atom line_to_atom(std::string line){
    std::istringstream iss(line);
    int elem;
    double x, y, z;
    iss >> elem >> x >> y >> z;
    Atom a;
    a.element = elem;
    a.location = arma::vec({x, y, z});
    a.beta = atom_lookup[elem][0];
    a.Z = atom_lookup[elem][1];
    return a;
}

// Function to read in atoms from xyz file path
std::vector<Atom> read_xyz_file(std::string atoms_file_path){
    // Define overall vector
    std::vector<Atom> atoms;

    // Open the file with fstream
    std::ifstream atoms_file(atoms_file_path);

    // Ignore both header lines
    std::string header_line_1;
    std::getline(atoms_file, header_line_1);
    std::string header_line_2;
    std::getline(atoms_file, header_line_2);

    // Iterate through the rest of the lines (while they exist)
    std::string line;
    while (std::getline(atoms_file, line)){
        // Get the line as an Atom struct
        Atom line_atom = line_to_atom(line);
        atoms.push_back(line_atom);
    }

    return atoms;
}

// Function to generate descending lexicographical triplets from L (overall angular momentum)
std::vector<std::vector<int>> generate_momentum_triplets(int L){
    std::vector<std::vector<int>> momentum_triplets;
    // Iterate downwards for the first element, from L to 0
    for (int i = L; i >= 0; i--){
        // Iterate downwards for the second element, from remaining sum down
        for (int j = (L-i); j >= 0; j--){
            // The third element is fixed based on the other two
            int k = L - i - j;
            // Define the triplet as a vector
            std::vector<int> triplet = {i, j, k};
            // Add to momentum triplets
            momentum_triplets.push_back(triplet);
        }
    }
    return momentum_triplets;
}

// Function to build basis functions from atoms
std::vector<Sto3G_Basis> build_basis_functions(const std::vector<Atom>& atoms){
    // Define a vector of basis functions to return
    std::vector<Sto3G_Basis> basis_functions;
    // Iterate over atoms
    for (int i = 0; i < atoms.size(); i++){
    //for (const auto& atom : atoms){
        Atom atom = atoms[i];
        // Get vector of shells for the atom
        std::vector<Sto3G_Shell> shells = shell_lookup.at(atom.element);
        // Iterate through shells
        for (const auto& shell : shells){
            // Generate angular momentum triplets for the given shell
            std::vector<std::vector<int>> shell_triplets = generate_momentum_triplets(shell.L);
            // Iterate through triplets
            for (const auto& triplet : shell_triplets){
                // Make a basis object to add to basis vector
                Sto3G_Basis basis = {atom.element, i, shell.ionization_energy, atom.location, triplet, shell};
                // Add this to overall basis function vector
                basis_functions.push_back(basis);
            }

        }

    }
    return basis_functions;
}

// Center of overlap
arma::vec center_of_overlap(arma::vec center_a, arma::vec center_b, double alpha, double beta){

    // Multiply center of a by alpha, center of b by beta, and add
    arma::vec numerator = (alpha * center_a) + (beta * center_b);
    // Divide this vector by alpha + beta to get center of product
    arma::vec center_p = numerator / (alpha + beta);

    return center_p;
}

// Define structure for 3D Gaussian Overlap
struct GaussOverlap3D{
    arma::vec center_a;
    arma::vec center_b;
    double alpha;
    double beta;
    std::vector<int> lA;
    std::vector<int> lB;
    // Calculate the center of the gaussian
    arma::vec center_p;
    GaussOverlap3D(arma::vec ca, arma::vec cb, double a, double b, std::vector<int> l1, std::vector<int> l2) : center_a(ca), center_b(cb), alpha(a), beta(b), lA(l1), lB(l2){
        center_p = center_of_overlap(ca, cb, a, b);
    }
};

// Helper function - evaluate factorial
int factorial(int n){
    // Error checking - just in case n < 0
    if (n < 0){
        throw std::invalid_argument("Input to factorial should never be negative");
    }
    // Set factorial to 1
    int result = 1;
    // Iterate down from n to 1
    for (int i = n; i > 1; i--){
        result *= i;
    }
    return result;
}

// Helper function - evaluate double factorial
int double_factorial(int n){
    // Return 1 if the integer is less than or equal to 0 (accounts for -1)
    if (n <= 0){
        return 1;
    }
    // Set the result to 1
    int result = 1;
    // Iterate through from n down to 2, with step of 2
    for (int i = n; i >= 2; i -= 2){
        result *= i;
    }
    return result;
}

// Helper function - evaluate the binomial
int evaluate_binomial(int upper, int lower){
    // Evaluate the numerator
    int numerator = factorial(upper);
    // Evaluate the denominator
    int denominator = factorial(lower) * factorial((upper - lower));
    return numerator / denominator;
}


// Calculate the prefactor to the full summation
double calculate_prefactor(const GaussOverlap3D& overlap, int comp){
    double diff = overlap.center_a[comp] - overlap.center_b[comp];
    double numerator = (overlap.alpha * overlap.beta) * std::pow(diff, 2.0);
    double denom = (overlap.alpha + overlap.beta);
    double term1 = std::exp((-numerator / denom));
    double term2 = std::sqrt((M_PI / denom));
    return term1 * term2;
}


// Perform the double summation loop for each x/y/z S_AB term
double perform_double_summation(const GaussOverlap3D& overlap, int comp){
    // Summation value to return
    double overall_sum = 0.0;
    // Overall summation, for index i
    for (int i = 0; i <= overlap.lA[comp]; i++){
        // Nested summation, for index j
        for (int j = 0; j <= overlap.lB[comp]; j++){
            // Check if i + j is even, and if so, include it (odd terms don't contribute)
            if ((i + j) % 2 == 0){
                // Calculate the ij double factorial term
                double ij_term = double_factorial((i + j - 1));
                // Calculate the a exponential term
                double a_exp = std::pow((overlap.center_p[comp] - overlap.center_a[comp]), (overlap.lA[comp] - i));
                // Calculate the b exponential term
                double b_exp = std::pow((overlap.center_p[comp] - overlap.center_b[comp]), (overlap.lB[comp] - j));
                // Calculate the denominator term
                double denom_term = std::pow((2.0 * (overlap.alpha + overlap.beta)), ((i + j) / 2.0));
                // Calculate multiplicative binomial prefactor
                double binom_prefactor = evaluate_binomial(overlap.lA[comp], i) * evaluate_binomial(overlap.lB[comp], j);
                // Combine it all together
                double summation_term = binom_prefactor * (ij_term * a_exp * b_exp) / denom_term;
                // Add to overall summation
                overall_sum += summation_term;
            }
        }
    }
    return overall_sum;
}

// Get the matrix S_AB term for a given gaussian overlap object
double get_overlap_matrix_term(const GaussOverlap3D& overlap){
    // Calculate X term
    double sab_x = calculate_prefactor(overlap, 0) * perform_double_summation(overlap, 0);
    // Calculate Y term
    double sab_y = calculate_prefactor(overlap, 1) * perform_double_summation(overlap, 1);
    // Calculate Z term
    double sab_z = calculate_prefactor(overlap, 2) * perform_double_summation(overlap, 2);

    // Return all three x/y/z multiplied together
    return (sab_x * sab_y * sab_z);
}

double get_normalization_term(const GaussOverlap3D& prim_self){
    // Calculate the overlap term for the primitive with itself (S_AA)
    double self_overlap = get_overlap_matrix_term(prim_self);
    // Return 1 / sqrt of S_AA
    double normalization_term = 1.0 / std::sqrt(self_overlap);
    return normalization_term;
}

// Function to get an S mu nu overlap term from two basis objects
double get_overlap_mu_nu_term(const Sto3G_Basis& b1, const Sto3G_Basis& b2){
    // We will sum over all primitives, start with 0
    double overlap_sum = 0.0;

    // Iterate through primitives of basis 1 (k)
    for (int k = 0; k < b1.sto3g.sto3g.size(); k++){
        // Iterate through primitives of basis 2 (l)
        for (int l = 0; l < b2.sto3g.sto3g.size(); l++){

            // Make a primitive 3D gaussian object for b1 with itself
            GaussOverlap3D prim_kk = GaussOverlap3D(b1.center, b1.center, b1.sto3g.sto3g[k].exp, b1.sto3g.sto3g[k].exp, b1.triplet, b1.triplet);
            // Make a primitive 3D gaussian object for b2 with itself
            GaussOverlap3D prim_ll = GaussOverlap3D(b2.center, b2.center, b2.sto3g.sto3g[l].exp, b2.sto3g.sto3g[l].exp, b2.triplet, b2.triplet);
            // Make a primitive 3D gaussian object for b1 with b2
            GaussOverlap3D prim_kl = GaussOverlap3D(b1.center, b2.center, b1.sto3g.sto3g[k].exp, b2.sto3g.sto3g[l].exp, b1.triplet, b2.triplet);

            // Get normalization terms from prim_kk and prim_ll
            double norm_kmu = get_normalization_term(prim_kk);
            double norm_lnu = get_normalization_term(prim_ll);
            // Get overlap between the two primitive shells (kl)
            double s_kl = get_overlap_matrix_term(prim_kl);
            // Get the contraction coefficients for this combo of k and l
            double d_kmu = b1.sto3g.sto3g[k].cc;
            double d_lnu = b2.sto3g.sto3g[l].cc;

            // Multiply all of these terms together to get the value for k, l
            double kl_sum_term = norm_kmu * norm_lnu * s_kl * d_kmu * d_lnu;
            // Add to overall summation
            overlap_sum += kl_sum_term;
        }
    }
    return overlap_sum;
}

// Function to populate the overlap matrix S
arma::mat populate_overlap_matrix(const std::vector<Sto3G_Basis>& basis_vec){
    // Make an empty arma::mat of size NxN
    arma::mat overlap_matrix(basis_vec.size(), basis_vec.size(), arma::fill::zeros);
    // Upper triangular iteration, first through entire vector
    for (int i = 0; i < basis_vec.size(); i++){
        // Iterate through remainder of vector starting at i (so diagonal is included)
        for (int j = i; j < basis_vec.size(); j++){
            // Get the term for that combination of bases
            double s_mu_nu = get_overlap_mu_nu_term(basis_vec[i], basis_vec[j]);
            // Store this in the matrix at both i,j and j,i (populates across diagonal)
            overlap_matrix(i, j) = s_mu_nu;
            overlap_matrix(j, i) = s_mu_nu;
        }
    }
    // Convert overlap matrix from AU to eV
    return overlap_matrix;
}

/*
CNDO/2 CODE
*/

// Function to compute the 0[0] term
double compute_zero_term(const arma::vec& rA, const arma::vec& rB, double sigmaA, double sigmaB){
    // Compute the UA term
    double UA = std::pow(M_PI * sigmaA, 1.5);
    // Compute the UB term
    double UB = std::pow(M_PI * sigmaB, 1.5);
    // Compute the V2 term
    double V2 = 1.0 / (sigmaA + sigmaB);
    // Compute distance between the two
    arma::vec r = rA - rB;
    double r2 = arma::dot(r, r);
    // Compute T term
    double T = r2 * V2;

    // If they are on top of each other (effectively r=t=0), ignore the erf part
    double zero_term;
    if (T < 1e-12){
        // Compute result minus error function
        zero_term = (UA * UB * 2.0) * std::sqrt((V2 / M_PI));
    }
    else{
        // Nominal case, compute result with error function
        double T_sqrt = std::sqrt(T);
        double multiplier = 0.5 * std::sqrt((M_PI / T)) * std::erf(T_sqrt);
        zero_term = UA * UB * multiplier * std::sqrt((2.0 * V2)) * std::sqrt((2.0 / M_PI));
    }
    return zero_term;
}

// Function to compute a gamma term (for a pair of atoms)
double compute_gamma(const Atom& atomA, const Atom& atomB){
    // Total gamma sum
    double gamma_sum = 0.0;
    // Get the s shell object for atom A
    Sto3G_Shell sshell_A = shell_lookup.at(atomA.element)[0];
    // Get the s shell object for atom B
    Sto3G_Shell sshell_B = shell_lookup.at(atomB.element)[0];
    // Four loops, over each of the shell's bases twice
    for (int k1 = 0; k1 < sshell_A.sto3g.size(); k1++){
        // Calculate normalization term for k1
        GaussOverlap3D prim_k1 = GaussOverlap3D(atomA.location, atomA.location, sshell_A.sto3g[k1].exp, sshell_A.sto3g[k1].exp, {0,0,0}, {0,0,0});
        double norm_k1 = get_normalization_term(prim_k1);

        for(int k2 = 0; k2 < sshell_A.sto3g.size(); k2++){
            // Calculate normalization term for k2
            GaussOverlap3D prim_k2 = GaussOverlap3D(atomA.location, atomA.location, sshell_A.sto3g[k2].exp, sshell_A.sto3g[k2].exp, {0,0,0}, {0,0,0});
            double norm_k2 = get_normalization_term(prim_k2);

            // Can calculate sigma A at this level, so we don't repeat it
            double sigma_A = 1.0 / (sshell_A.sto3g[k1].exp + sshell_A.sto3g[k2].exp);

            for (int l1 = 0; l1 < sshell_B.sto3g.size(); l1++){
                // Calculate normalization term for l1
                GaussOverlap3D prim_l1 = GaussOverlap3D(atomB.location, atomB.location, sshell_B.sto3g[l1].exp, sshell_B.sto3g[l1].exp, {0,0,0}, {0,0,0});
                double norm_l1 = get_normalization_term(prim_l1);

                for (int l2 = 0; l2 < sshell_B.sto3g.size(); l2++){
                    // Calculate normalization term for l2
                    GaussOverlap3D prim_l2 = GaussOverlap3D(atomB.location, atomB.location, sshell_B.sto3g[l2].exp, sshell_B.sto3g[l2].exp, {0,0,0}, {0,0,0});
                    double norm_l2 = get_normalization_term(prim_l2);

                    // Calculate sigma B
                    double sigma_B = 1.0 / (sshell_B.sto3g[l1].exp + sshell_B.sto3g[l2].exp);

                    // Get the zero term for this combo
                    double current_zero_term = compute_zero_term(atomA.location, atomB.location, sigma_A, sigma_B);

                    // Calculate total normalization term
                    double norm_total = norm_k1 * norm_k2 * norm_l1 * norm_l2;

                    // Add to total gamma sum with contraction coefficients
                    gamma_sum += (norm_total * sshell_A.sto3g[k1].cc * sshell_A.sto3g[k2].cc * sshell_B.sto3g[l1].cc * sshell_B.sto3g[l2].cc * current_zero_term);
                }
            }
        }
    }
    return gamma_sum;
}

// Function to compute full gamma matrix (over all atoms)
arma::mat compute_gamma_matrix(const std::vector<Atom>& atoms){
    arma::mat gamma_matrix(atoms.size(), atoms.size(), arma::fill::zeros);
    // Iterate over atoms, as upper diagonal
    for (int i = 0; i < atoms.size(); i++){
        for(int j = i; j < atoms.size(); j++){
            // Get the atom A, B for this location
            Atom atomA = atoms[i];
            Atom atomB = atoms[j];
            double gamma_ij = compute_gamma(atomA, atomB);
            // Populate gamma matrix across diagonals (will duplicate diagonal)
            gamma_matrix(i,j) = gamma_ij;
            gamma_matrix(j,i) = gamma_ij;
        }
    }
    // Convert gamma to eV (from AU)
    arma::mat gamma_matrix_ev = gamma_matrix * 27.211324570273;
    return gamma_matrix_ev;
}

// Helper function, compute off diagonal term of fock matrix
double compute_fock_off_diag_term(double betaA, double betaB, double s_mu_nu, double p_mu_nu, double gamma_AB, bool full_fock = true){
    // Calculate term 1
    double term1 = 0.5 * (betaA + betaB) * s_mu_nu;
    // Calculate total (depending on full vs core fock)
    double result;
    if (full_fock){
        // Calculate term 2 if full fock
        double term2 = p_mu_nu * gamma_AB;
        result = (term1 - term2);
    }
    else{
        // If just core hamiltonian, don't include second portion
        result = term1;
    }

    return result;
}

// Helper function, get total density for a given atom
double compute_per_atom_P(const arma::mat& P, const std::vector<Sto3G_Basis>& basis_vec, int atom_index){
    // Total sum of density
    double P_sum = 0.0;
    // Iterate over basis functions
    for (int i = 0; i < basis_vec.size(); i++){
        // Get the basis object
        Sto3G_Basis basis = basis_vec[i];
        // If the atom corresponding to basis is atomic index
        if (basis.atomic_index == atom_index){
            // Add the density for that atom to sum
            P_sum += P(i,i); // we sum the diagonal elements
        }
    }
    return P_sum;
}

// Helper function, compute diagonal term of fock matrix
double compute_fock_diag_term(const std::vector<Sto3G_Basis>& bases, const std::vector<Atom>& atoms, const arma::mat& gamma_matrix, const arma::mat& P_tot, int atom_index, double IA_term, double p_nu_nu, bool full_fock = true){
    // Sum term gamma AC, total density on all atoms except A
    double gamma_AC_sum = 0.0;
    // Gamma AA, gamma term for atom A
    double gamma_AA;
    // Density term PAA, total density on atom A
    double P_tot_AA;
    // Z value for atom A
    double Z_A;

    // Iterate through atoms to simultaneously determine individual atom and sum parameters
    for (int i = 0; i < atoms.size(); i++){
        // If the index is atomic index, it's the "atom in question"
        if (i == atom_index){
            if (full_fock){
                // Calculate total density for atom A (if full fock)
                P_tot_AA = compute_per_atom_P(P_tot, bases, i);
            }
            // If atom in question, set variables for its gamma and Z
            gamma_AA = gamma_matrix(i,i);
            Z_A = atoms[i].Z;
        }
        // Otherwise, add the other atomic parameters to sum term
        else{
            if (full_fock){
                // Calculate the density sum for atom CC
                double P_tot_CC = compute_per_atom_P(P_tot, bases, i);
                // Add to the AC summation term, along with Z and gamma
                gamma_AC_sum += ((P_tot_CC - atoms[i].Z) * gamma_matrix(atom_index, i));
            }
            else{
                // Case where we are calculating core hamiltonian, no total density
                gamma_AC_sum -= (atoms[i].Z * gamma_matrix(atom_index, i));
            }
        }
    }
    // Compute the second term in full equation (different for core vs full hamiltonian)
    double second_term;
    if (full_fock){
        second_term = ((P_tot_AA - Z_A) - (p_nu_nu - 0.5)) * gamma_AA;
    }
    else{
        second_term = - (Z_A - 0.5) * gamma_AA;
    }

    // Calculate total full fock nu nu term and return
    double fock_nu_nu = -IA_term + second_term + gamma_AC_sum;

    return fock_nu_nu;
}

// Function to compute fock matrix
arma::mat compute_fock_matrix(const std::vector<Sto3G_Basis>& basis_vec, const std::vector<Atom>& atoms, const arma::mat& P_ab, const arma::mat& P_tot, const arma::mat& gamma_matrix, const arma::mat& S, bool full_fock = true, arma::vec elec_field = arma::zeros<arma::vec>(3)){
    // Make an empty arma::mat of size NxN
    arma::mat fock_matrix(basis_vec.size(), basis_vec.size(), arma::fill::zeros);

    // Iterate over basis vector as upper diagonal
    for (int mu = 0; mu < basis_vec.size(); mu++){
        for (int nu = mu; nu < basis_vec.size(); nu++){
            // Get each of the bases
            Sto3G_Basis basis_mu = basis_vec[mu];
            Sto3G_Basis basis_nu = basis_vec[nu];
            // Get atoms corresponding to bases
            Atom atom_A = atoms[basis_mu.atomic_index];
            Atom atom_B = atoms[basis_nu.atomic_index];
            // Case 1 - off diagonal term
            if (mu != nu){
                // Extract the gamma AB term
                double gamma_AB = gamma_matrix(basis_mu.atomic_index, basis_nu.atomic_index);
                // Compute off diagonal fock term
                double fock_mu_nu = compute_fock_off_diag_term(atom_A.beta, atom_B.beta, S(mu, nu), P_ab(mu, nu), gamma_AB, full_fock);
                // Add to fock matrix (symmetric)
                fock_matrix(mu, nu) = fock_mu_nu;
                fock_matrix(nu, mu) = fock_mu_nu;
            }
            // Case 2 - diagonal term
            else{
                // Compute diagonal fock term
                double fock_nu_nu = compute_fock_diag_term(basis_vec, atoms, gamma_matrix, P_tot, basis_nu.atomic_index, basis_nu.sto3g.IA_coeff, P_ab(nu, nu), full_fock);
                // Compute field shift for diagonal term (zero nominally)
                double electric_field_shift = arma::dot(elec_field, atom_A.location);
                // Add this to the diagonal term
                fock_nu_nu += electric_field_shift;
                // Add to fock matrix
                fock_matrix(nu, nu) = fock_nu_nu;
            }
        }
    }
    // Return matrix
    return fock_matrix;
}

// Structure for eigenvalue solution
struct EigenSolutionStandard{
    arma::mat C;
    arma::vec E;
};

// Structure for SCF solution
struct SCFSolution{
    arma::mat F_alpha;
    arma::mat F_beta;
    arma::mat P_alpha;
    arma::mat P_beta;
    arma::mat H;
    arma::vec E_alpha;
    arma::vec E_beta;
    arma::mat C_alpha;
    arma::mat C_beta;
};

// Helper function, to solve standard eigenvalue problem
EigenSolutionStandard solve_standard_eigenvalue(arma::mat F){
    // Define variables to solve
    arma::vec E; // Energies
    arma::mat C; // Coefficients

     // Solve the standard eigenvalue problem F C = C E
    bool eigval_success = arma::eig_sym(E, C, F);
    if (! eigval_success){
        throw std::runtime_error("Error, eigen decomposition of fock matrix failed.");
    }

    EigenSolutionStandard eigsol;
    eigsol.E = E;
    eigsol.C = C;

    return eigsol;
}

// Function to calculate new P matrices from C
arma::mat calculate_P_from_C(const arma::mat& C, int n_elec){
    // Initialize empty P matrix of dimensions of C
    arma::mat P(C.n_rows, C.n_cols);
    // Loop over mu dimension of C
    for (int mu = 0; mu < C.n_rows; mu++){
        // Loop over nu dimension of C
        for (int nu = 0; nu < C.n_cols; nu++){
            // Loop over the number of electrons, sum it up
            double mu_nu_sum = 0.0;
            for (int i = 0; i < n_elec; i++){
                // Get the mu nu term for this electron, add to sum
                mu_nu_sum += C(mu, i) * C(nu, i);
            }
            // Place this element in the P matrix at mu, nu
            P(mu, nu) = mu_nu_sum;
        }
    }
    return P;
}

// Overall function to implement SCF
SCFSolution cndo2_scf(const std::vector<Sto3G_Basis>& basis_vec, const std::vector<Atom>& atoms, const arma::mat& gamma_matrix, const arma::mat& S, int p, int q, arma::vec ext_field = arma::zeros<arma::vec>(3), double tolerance = 1e-6, int max_iterations = 1000){
    // Initial guess for P matrices, zeros, dimension of bases
    arma::mat P_alpha(basis_vec.size(), basis_vec.size(), arma::fill::zeros);
    arma::mat P_beta(basis_vec.size(), basis_vec.size(), arma::fill::zeros);
    arma::mat P_tot = P_alpha + P_beta;

    // Variables to set in loop
    arma::mat fock_alpha, fock_beta;
    arma::mat P_alpha_old, P_beta_old, P_tot_old;
    arma::mat C_alpha, C_beta;
    arma::vec E_alpha, E_beta;

    // Loop conditionals
    bool converged = false;
    int iteration_number = 0;

    // Full loop, while not converged
    while (!converged){

        // Build the fock matrices
        fock_alpha = compute_fock_matrix(basis_vec, atoms, P_alpha, P_tot, gamma_matrix, S, true, ext_field);
        fock_beta = compute_fock_matrix(basis_vec, atoms, P_beta, P_tot, gamma_matrix, S, true, ext_field);

        // Solve the eigenvalue problem to get C's, E's
        EigenSolutionStandard fock_sol_alpha = solve_standard_eigenvalue(fock_alpha);
        EigenSolutionStandard fock_sol_beta = solve_standard_eigenvalue(fock_beta);

        // Extract from solution object
        C_alpha = fock_sol_alpha.C;
        C_beta = fock_sol_beta.C;
        E_alpha = fock_sol_alpha.E;
        E_beta = fock_sol_beta.E;

        // Copy old P matrices over
        P_alpha_old = P_alpha;
        P_beta_old = P_beta;
        P_tot_old = P_tot;

        // Use C's to construct new P matrices
        P_alpha = calculate_P_from_C(C_alpha, p);
        P_beta = calculate_P_from_C(C_beta, q);
        P_tot = P_alpha + P_beta;

        // Check convergence of each
        bool alpha_converged = arma::approx_equal(P_alpha, P_alpha_old, "absdiff", tolerance);
        bool beta_converged = arma::approx_equal(P_beta, P_beta_old, "absdiff", tolerance); 

        // If both converged, we end loop
        if (alpha_converged && beta_converged){
            converged = true;
        }

        // Also check if it maxed out on iterations
        if (iteration_number > max_iterations){
            throw std::runtime_error("SCF loop failed to converge, exceeded maximum iterations.");
        }
        // Update the iteration number
        iteration_number++;
    }

    // Construct the custom solution object to return
    SCFSolution solution;
    solution.F_alpha = compute_fock_matrix(basis_vec, atoms, P_alpha, P_tot, gamma_matrix, S, true, ext_field);
    solution.F_beta = compute_fock_matrix(basis_vec, atoms, P_beta, P_tot, gamma_matrix, S, true, ext_field);
    solution.P_alpha = P_alpha;
    solution.P_beta = P_beta;
    solution.H = compute_fock_matrix(basis_vec, atoms, P_alpha, P_tot, gamma_matrix, S, false, ext_field);
    solution.C_alpha = C_alpha;
    solution.E_alpha = E_alpha;
    solution.C_beta = C_beta;
    solution.E_beta = E_beta;

    return solution;

}

// Structure to store and return all energy components
struct EnergySolution{
    double E_electronic;
    double E_nuclear;
    double E_total;
};

// Function to calculate total fock energy
EnergySolution calculate_fock_energy(const SCFSolution& fock_solution, const std::vector<Atom>& atoms){
    // Calculate the nuclear repulsion energy
    double nuclear_repulsion_energy = 0.0;
    // Iterate over atom pairs
    for (int i = 0; i < atoms.size(); i++){
        for (int j = i+1; j < atoms.size(); j++){
            // Atom pair
            Atom atom_A = atoms[i];
            Atom atom_B = atoms[j];
            // Get the radius between the atoms
            arma::vec r = atom_A.location - atom_B.location;
            double r_AB = arma::norm(r);
            // Add current nuclear repulsion to total
            nuclear_repulsion_energy += ((atom_A.Z * atom_B.Z) / r_AB);
        }
    }
    // Convert nuclear repulsion energy to eV (from AU)
    nuclear_repulsion_energy *= 27.211324570273;

    // Calculate the alpha energy term
    arma::mat total_alpha_mat = fock_solution.P_alpha % (fock_solution.H + fock_solution.F_alpha);
    double alpha_energy_term = 0.5 * arma::sum(arma::sum(total_alpha_mat));

    // Calculate the beta energy term
    arma::mat total_beta_mat = fock_solution.P_beta % (fock_solution.H + fock_solution.F_beta);
    double beta_energy_term = 0.5 * arma::sum(arma::sum(total_beta_mat));

    // Total sum for total energy
    double fock_energy = alpha_energy_term + beta_energy_term + nuclear_repulsion_energy;

    // Define energy structure to return
    EnergySolution esol;
    esol.E_electronic = alpha_energy_term + beta_energy_term;
    esol.E_nuclear = nuclear_repulsion_energy;
    esol.E_total = fock_energy;

    return esol;
}