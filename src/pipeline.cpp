#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <iomanip>

namespace fs = std::filesystem;
using json = nlohmann::json;

#include "polarity.hpp"
#include "cndo2_utils.hpp"

// Structure to save results, to ultimately sort
struct PropertyResult {
    std::string molecule;
    double value;
};

// Helper function - get molecule name from path (i.e. atoms/CH4.xyz -> CH4)
std::string molecule_name_from_path(const fs::path& atoms_file_path) {
    return atoms_file_path.stem().string();
}

// Write out a csv file, sorting by values in struct
// Note - highest to lowest sorting!
void write_sorted_csv(
    const std::string& output_path,
    const std::string& value_column_name,
    std::vector<PropertyResult> results
) {
    std::sort(
        results.begin(),
        results.end(),
        [](const PropertyResult& a, const PropertyResult& b) {
            return a.value > b.value;
        }
    );

    std::ofstream out(output_path);
    if (!out) {
        throw std::runtime_error("Could not open output file: " + output_path);
    }

    out << "molecule," << value_column_name << "\n";
    out << std::setprecision(12);

    for (const auto& result : results) {
        out << result.molecule << "," << result.value << "\n";
    }
}

int main() {
    std::string input_path = "./input";
    std::string output_dir = "./pipeline_outputs";

    // Create the output directory if it does not exist
    fs::create_directories(output_dir);

    // Vectors to store the results in, to ultimately sort
    std::vector<PropertyResult> dipole_original;
    std::vector<PropertyResult> iso_original;
    std::vector<PropertyResult> aniso_original;

    std::vector<PropertyResult> dipole_ao;
    std::vector<PropertyResult> iso_ao;
    std::vector<PropertyResult> aniso_ao;

    for (const auto& entry : fs::directory_iterator(input_path)) {
        if (entry.path().extension() != ".json") {
            continue;
        }

        std::ifstream config_file(entry.path());
        if (!config_file) {
            throw std::runtime_error("Could not open config file: " + entry.path().string());
        }

        json config = json::parse(config_file);

        fs::path atoms_file_path = config["atoms_file_path"].get<std::string>();
        int p = config["num_alpha_electrons"];
        int q = config["num_beta_electrons"];

        // Extract molecule name
        std::string molecule = molecule_name_from_path(atoms_file_path);

        // Print molecule
        std::cout << "Computing results for " << molecule << std::endl;

        // Original atom centered method
        arma::vec dipole_vector = compute_dipole_from_xyz(atoms_file_path.string(), p, q);
        double dipole = dipole_mag(dipole_vector);

        arma::mat alpha = calculate_polarizability_tensor(atoms_file_path.string(), p, q, false);
        double iso_polarizability = isotropic_polarizability(alpha);
        double aniso_polarizability = anisotropic_polarizability(alpha);

        // Dipole AO integral method
        arma::vec dipole_vector_ao = compute_ao_dipole_from_xyz(atoms_file_path.string(), p, q);
        double dipole_new = dipole_mag(dipole_vector_ao);

        arma::mat alpha_ao_mat = calculate_polarizability_tensor(atoms_file_path.string(), p, q, true);
        double iso_polarizability_new = isotropic_polarizability(alpha_ao_mat);
        double aniso_polarizability_new = anisotropic_polarizability(alpha_ao_mat);

        dipole_original.push_back({molecule, dipole});
        iso_original.push_back({molecule, iso_polarizability});
        aniso_original.push_back({molecule, aniso_polarizability});

        dipole_ao.push_back({molecule, dipole_new});
        iso_ao.push_back({molecule, iso_polarizability_new});
        aniso_ao.push_back({molecule, aniso_polarizability_new});
    }

    // Write all of the results out to csv files, original and updated method
    write_sorted_csv(output_dir + "/dipole_original.csv", "dipole_original", dipole_original);
    write_sorted_csv(output_dir + "/iso_polarizability_original.csv", "iso_polarizability_original", iso_original);
    write_sorted_csv(output_dir + "/aniso_polarizability_original.csv", "aniso_polarizability_original", aniso_original);
    write_sorted_csv(output_dir + "/dipole_ao_integral.csv", "dipole_ao_integral", dipole_ao);
    write_sorted_csv(output_dir + "/iso_polarizability_ao_integral.csv", "iso_polarizability_ao_integral", iso_ao);
    write_sorted_csv(output_dir + "/aniso_polarizability_ao_integral.csv", "aniso_polarizability_ao_integral", aniso_ao);

    // Print that pipeline finished
    std::cout << "Pipeline complete. Results written to " << output_dir << "\n";

    return 0;
}