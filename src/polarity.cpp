#include <armadillo>
#include <cmath>
#include "cndo2_utils.hpp"

/*
PERMANENT DIPOLE CALCULATION
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
arma::vec compute_dipole_from_xyz(std::string atoms_file_path, int p, int q){

    // Read atoms from the file
    std::vector<Atom> atoms = read_xyz_file(atoms_file_path);

    // Get vector of basis functions
    std::vector<Sto3G_Basis> bases = build_basis_functions(atoms);

    // Get the overlap matrix
    arma::mat S = populate_overlap_matrix(bases);

    // Get the gamma matrix
    arma::mat gamma_mat = compute_gamma_matrix(atoms);

    // Run SCF
    SCFSolution scf_sol = cndo2_scf(bases, atoms, gamma_mat, S, p, q);

    // Calculate dipole from SCF solution
    arma::vec dipole = calculate_permanent_dipole(atoms, bases, scf_sol);

    return dipole;
}

// Helper - norm of dipole
double dipole_norm(const arma::vec& dipole_vec){
    double dipole = arma::norm(dipole_vec);
    return dipole;
}

/*
PROPERTY CALCULATION
*/

//assuming dipole moment is arma::vec of length 3 = mu
//assuming polarizability tensor is (3 x 3) arma::mat = alpha
double dipole_mag(const arma::vec& mu) {
    return std::sqrt(arma::dot(mu, mu));
}

double iso_polarizability(const arma::mat& alpha){
    return (alpha(0,0) + alpha(1,1) + alpha(2,2)) / 3.0;
}

double aniso_polarizability(const arma::mat& alpha){
    double xy_delta = (alpha(0,0) - alpha(1,1)) * (alpha(0,0) - alpha(1,1));
    double yz_delta = (alpha(1,1) - alpha(2,2)) * (alpha(1,1) - alpha(2,2));
    double xz_delta = (alpha(0,0) - alpha(2,2)) * (alpha(0,0) - alpha(2,2));

    return std::sqrt(0.5 * (xy_delta + yz_delta + xz_delta));
}