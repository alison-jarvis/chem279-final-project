#pragma once

#include <armadillo>
#include <cmath>
#include "cndo2_utils.hpp"

/*
PERMANENT DIPOLE CALCULATION
*/

// Net electron population
std::vector<double> calculate_NA(const std::vector<Atom>& atoms, const std::vector<Sto3G_Basis>& bases, const arma::mat& P_alpha, const arma::mat& P_beta);

// Net atomic charge
std::vector<double> calculate_qA(const std::vector<Atom>& atoms, const std::vector<double>& populations);

// Compute dipole
arma::vec calculate_permanent_dipole(const std::vector<Atom>& atoms, const std::vector<Sto3G_Basis>& bases, const SCFSolution& scf_sol);

// Wrapper, dipole from config options
arma::vec compute_dipole_from_xyz(
    std::string atoms_file_path, 
    int p, 
    int q, 
    arma::vec ext_electric_field = arma::zeros<arma::vec>(3)
);

/*
POLARIZABILITY CALCULATION
*/

arma::mat calculate_polarizability_tensor(std::string atoms_file_path, int p, int q, bool with_integral = true, double f = 1.0e-4);


/*
PROPERTY CALCULATION
*/

double dipole_mag(const arma::vec& mu);

double isotropic_polarizability(const arma::mat& alpha);

double anisotropic_polarizability(const arma::mat& alpha);