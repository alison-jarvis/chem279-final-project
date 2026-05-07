#pragma once

#include <armadillo>
#include <vector>
#include <string>
#include <map>
#include <nlohmann/json.hpp>

/*
Data structures
*/

struct Atom{
    int element;
    arma::vec location;
    double beta;
    int Z;
};

struct Primitive_Gauss{
    double exp;
    double cc;
};

struct Sto3G_Shell{
    int L;
    double ionization_energy;
    double IA_coeff;
    std::vector<Primitive_Gauss> sto3g;

    Sto3G_Shell(int L, double ie, double iac, std::string sto3g_path);
};

struct Sto3G_Basis{
    int element;
    int atomic_index;
    double ionization_energy;
    arma::vec center;
    std::vector<int> triplet;
    Sto3G_Shell sto3g;
};

struct EigenSolutionStandard{
    arma::mat C;
    arma::vec E;
};

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


std::vector<Atom> read_xyz_file(std::string atoms_file_path);


std::vector<Sto3G_Basis> build_basis_functions(const std::vector<Atom>& atoms);


arma::mat populate_overlap_matrix(const std::vector<Sto3G_Basis>& basis_vec);

/*
PERMANENT DIPOLE AO INTEGRAL CALCULATION
*/

// Helper - get primitive dipole component
double get_primitive_dipole_component(const Sto3G_Basis& b1, const Sto3G_Basis& b2, int k, int prim1, int prim2);

// Helper get dipole component for two full AOs
double get_ao_dipole_component(const Sto3G_Basis& b1, const Sto3G_Basis& b2, int k);

// Helper - build the dipole matrix for given comp k, indexed by orbitals
arma::mat build_ao_dipole_matrix(const std::vector<Sto3G_Basis>& bases, int k);

// Wrapper - calculate permanent dipole from atoms, bases, SCF solution
arma::vec calculate_permanent_dipole_ao_integrals(const std::vector<Atom>& atoms, const std::vector<Sto3G_Basis>& bases, const SCFSolution& scf_sol);

// Wrapper - from xyz compute AO form of dipole
arma::vec compute_ao_dipole_from_xyz(std::string atoms_file_path, int p, int q, arma::vec ext_electric_field = arma::zeros<arma::vec>(3));

/*
CNDO/2 Things
*/
double compute_per_atom_P(const arma::mat& P, const std::vector<Sto3G_Basis>& basis_vec, int atom_index);

arma::mat compute_gamma_matrix(const std::vector<Atom>& atoms);

arma::mat compute_fock_matrix(
    const std::vector<Sto3G_Basis>& basis_vec, 
    const std::vector<Atom>& atoms, 
    const arma::mat& P_ab, 
    const arma::mat& P_tot, 
    const arma::mat& gamma_matrix, 
    const arma::mat& S, 
    bool full_fock = true, 
    arma::vec elec_field = arma::zeros<arma::vec>(3)
);

SCFSolution cndo2_scf(
    const std::vector<Sto3G_Basis>& basis_vec, 
    const std::vector<Atom>& atoms, 
    const arma::mat& gamma_matrix, 
    const arma::mat& S, 
    int p, 
    int q, 
    arma::vec ext_field = arma::zeros<arma::vec>(3), 
    double tolerance = 1e-6, 
    int max_iterations = 1000
);

struct EnergySolution{
    double E_electronic;
    double E_nuclear;
    double E_total;
};

EnergySolution calculate_fock_energy(
    const SCFSolution& fock_solution,
    const std::vector<Atom>& atoms
);