#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <filesystem>

namespace fs = std::filesystem;
using json = nlohmann::json;

#include "polarity.hpp"

int main() {

    std::string path = "./input";
    for (const auto& entry : fs::directory_iterator(path)) {
        // Parse config file path into json config
        std::ifstream config_file(entry.path());
        json config = json::parse(config_file);

        // Extract important information from config
        fs::path atoms_file_path = config["atoms_file_path"];
        int p = config["num_alpha_electrons"];
        int q = config["num_beta_electrons"];
        
        // Calculate dipole from molecule of interest
        arma::vec dipole_vector = compute_dipole_from_xyz(atoms_file_path, p, q);
        double dipole = dipole_mag(dipole_vector);

        // Calculate polarizability tensor from molecule of interest
        arma::mat alpha = calculate_polarizability_tensor(atoms_file_path, p, q);

        // Calculate derived properties from this
        double iso_polarizability = isotropic_polarizability(alpha);
        double aniso_polarizability = anisotropic_polarizability(alpha);

        std::cout << "File: " << atoms_file_path << std::endl;
        std::cout << "Permanent dipole magnitude: " << dipole << std::endl;
        std::cout << "Iso polarizability: " << iso_polarizability << std::endl;
        std::cout << "Aniso polarizability: " << aniso_polarizability << std::endl << std::endl;
    }

    return 0;
}
