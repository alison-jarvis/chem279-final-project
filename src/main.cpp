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
        double dipole = dipole_norm(dipole_vector);

        std::cout << "Permanent Dipole for " << atoms_file_path << ":" << std::endl << dipole_vector << std::endl;
        std::cout << "Magnitude: " << dipole << std::endl;
    }

    return 0;
}
