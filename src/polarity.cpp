#include <armadillo>
#include <cmath>
#include "cndo2_utils.hpp"

/*
ORIGINAL PERMANENT DIPOLE CALCULATION
*/

// Helper function - net electron population on atom A
std::vector<double> calculate_NA(const std::vector<Atom>& atoms, const std::vector<Sto3G_Basis>& bases, const arma::mat& P_alpha, const arma::mat& P_beta){
    // Calculate total P_mu nu
    arma::mat P_tot = P_alpha + P_beta;
    // Initialize vector of size atoms to store NA's
    std::vector<double> populations;
    // Iterate through atoms
    for (int i = 0; i < atoms.size(); i++){
        // Calculate density on that atom
        double P_tot_AA = compute_per_atom_P(P_tot, bases, i);
        // Add density to vector
        populations.push_back(P_tot_AA);
    }
    return populations;
}

// Helper function - net atomic charge
std::vector<double> calculate_qA(const std::vector<Atom>& atoms, const std::vector<double>& populations){
    // Initialize empty vector
    std::vector<double> atomic_charges;
    // Iterate through atoms
    for (int i = 0; i < atoms.size(); i++){
        // Compute q_A as Z_A - N_A
        double qA = static_cast<double>(atoms[i].Z) - populations[i];
        // Add to vector
        atomic_charges.push_back(qA);
    }
    return atomic_charges;
}

// Net function - compute dipole
arma::vec calculate_permanent_dipole(const std::vector<Atom>& atoms, const std::vector<Sto3G_Basis>& bases, const SCFSolution& scf_sol){
    // Calculate net electron populations
    std::vector<double> populations = calculate_NA(atoms, bases, scf_sol.P_alpha, scf_sol.P_beta);
    // Calculate net atomic charges
    std::vector<double> atomic_charges = calculate_qA(atoms, populations);
    // Initialize dipoles
    arma::vec dipoles = arma::zeros<arma::vec>(3);
    // Iterate through atoms and compute dipole
    for (int i = 0; i < atoms.size(); i++){
        // Add elements to dipole
        dipoles = dipoles + (atoms[i].location * atomic_charges[i]);
    }
    return dipoles;
}

// Wrapper - from xyz file, run SCF, compute dipole
arma::vec compute_dipole_from_xyz(std::string atoms_file_path, int p, int q, arma::vec ext_electric_field = arma::zeros<arma::vec>(3)){

    // Read atoms from the file
    std::vector<Atom> atoms = read_xyz_file(atoms_file_path);

    // Get vector of basis functions
    std::vector<Sto3G_Basis> bases = build_basis_functions(atoms);

    // Get the overlap matrix
    arma::mat S = populate_overlap_matrix(bases);

    // Get the gamma matrix
    arma::mat gamma_mat = compute_gamma_matrix(atoms);

    // Run SCF
    SCFSolution scf_sol = cndo2_scf(bases, atoms, gamma_mat, S, p, q, ext_electric_field);

    // Calculate dipole from SCF solution
    arma::vec dipole = calculate_permanent_dipole(atoms, bases, scf_sol);

    return dipole;
}

/*
POLARIZABILITY CALCULATION
*/

// Function to calculate polarizability tensor
arma::mat calculate_polarizability_tensor(std::string atoms_file_path, int p, int q, bool with_integral = true, double f = 1.0e-4){
    // Initialize tensor
    arma::mat alpha(3, 3, arma::fill::zeros);

    // Loop over field direction, j
    for (int j = 0; j < 3; j++){
        // Get positive external field in ith component
        arma::vec pos_ext_field = arma::zeros<arma::vec>(3);
        pos_ext_field(j) = f;
        // Calculate dipole for positive f
        arma::vec dipole_positive_f;
        if (with_integral){
            dipole_positive_f = compute_ao_dipole_from_xyz(atoms_file_path, p, q, pos_ext_field);
        }
        else{
            dipole_positive_f = compute_dipole_from_xyz(atoms_file_path, p, q, pos_ext_field);
        }

        // Get negative external field in ith component
        arma::vec neg_ext_field = arma::zeros<arma::vec>(3);
        neg_ext_field(j) = -f;
        // Calculate dipole for negative f
        arma::vec dipole_negative_f;
        if (with_integral){
            dipole_negative_f = compute_ao_dipole_from_xyz(atoms_file_path, p, q, neg_ext_field);
        }
        else{
            dipole_negative_f = compute_dipole_from_xyz(atoms_file_path, p, q, neg_ext_field);
        }

        // Loop over dipole components, i
        for (int i = 0; i < 3; i++){
            // Calculate the component with central differences
            double alpha_ij = (dipole_positive_f(i) - dipole_negative_f(i))/(2.0 * f);
            // Store in tensor
            alpha(i, j) = alpha_ij;
        }
    }
    return alpha;
}

/*
PROPERTY CALCULATION
*/

//assuming dipole moment is arma::vec of length 3 = mu
//assuming polarizability tensor is (3 x 3) arma::mat = alpha
double dipole_mag(const arma::vec& mu) {
    return std::sqrt(arma::dot(mu, mu));
}

// Isotropic polarizability
double isotropic_polarizability(const arma::mat& alpha){
    return (alpha(0,0) + alpha(1,1) + alpha(2,2)) / 3.0;
}

// Ansiotropic polarizability
double ansiotropic_polarizability(const arma::mat& alpha){

    // Symmetrize the off diagonal terms
    arma::mat a = 0.5 * (alpha + alpha.t());

    // Numerator diagonal terms
    double xy_delta = (a(0,0) - a(1,1)) * (a(0,0) - a(1,1));
    double yz_delta = (a(1,1) - a(2,2)) * (a(1,1) - a(2,2));
    double xz_delta = (a(0,0) - a(2,2)) * (a(0,0) - a(2,2));
    double numerator_term1 = xy_delta + yz_delta + xz_delta;

    // Numerator off diagonal terms
    double cross_term_sum = (a(0,1)*a(0,1)) + (a(1,2)*a(1,2)) + (a(0,2)*a(0,2));

    // Total term in brackets
    double combined = (numerator_term1 + (6.0 * cross_term_sum)) / 2.0;

    return std::sqrt(combined);
}