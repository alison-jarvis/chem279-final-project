#include <armadillo>
#include <cmath>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <filesystem>
#include <vector>
#include <map>
#include <iomanip>

namespace fs = std::filesystem;
using json = nlohmann::json;

#include "polarity.hpp"
#include "cndo2_utils.hpp"

// Helper function - get molecule name from path (i.e. atoms/CH4.xyz -> CH4)
std::string molecule_name_from_path(const fs::path& atoms_file_path) {
    return atoms_file_path.stem().string();
}


// Helper function - write values to csv file
// Field magnitude in first column, other columns molecule name, values
void write_field_dependence_csv(const std::string& output_path, 
    const std::vector<double>& field_values,
    const std::vector<std::string>& molecule_names,
    const std::map<std::string, std::vector<double>>& molecule_to_values
) {
    std::ofstream out(output_path);

    if (!out) {
        throw std::runtime_error("Could not open output file: " + output_path);
    }

    out << std::setprecision(12);

    // Header
    out << "field_magnitude";

    for (const auto& molecule : molecule_names) {
        out << "," << molecule;
    }

    out << "\n";

    // Rows
    for (size_t i = 0; i < field_values.size(); ++i) {

        out << field_values[i];

        for (const auto& molecule : molecule_names) {
            out << "," << molecule_to_values.at(molecule)[i];
        }

        out << "\n";
    }
}

int main() {

    std::string input_dir = "./input";
    std::string output_dir = "./testing_outputs";

    // Create the outputs directory if it doesn't exist
    fs::create_directories(output_dir);

    // Field values, in reasonable log magnitude range
    std::vector<double> field_values = {1.0e-6, 3.0e-6, 1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3};

    // Molecule order for csv columns
    std::vector<std::string> molecule_names;

    // Storage, maps molecule to vector of values (over field strength)
    std::map<std::string, std::vector<double>> iso_original;
    std::map<std::string, std::vector<double>> aniso_original;

    std::map<std::string, std::vector<double>> iso_ao;
    std::map<std::string, std::vector<double>> aniso_ao;

    // Loop through molecules
    for (const auto& entry : fs::directory_iterator(input_dir)) {

        if (entry.path().extension() != ".json") {
            continue;
        }

        std::ifstream config_file(entry.path());

        if (!config_file) {
            throw std::runtime_error(
                "Could not open config file: " +
                entry.path().string()
            );
        }

        json config = json::parse(config_file);

        fs::path atoms_file_path =
            config["atoms_file_path"].get<std::string>();

        int p = config["num_alpha_electrons"];
        int q = config["num_beta_electrons"];

        // Get molecule name and save in names vector
        std::string molecule = molecule_name_from_path(atoms_file_path);
        molecule_names.push_back(molecule);

        // Initialize vectors
        iso_original[molecule] = {};
        aniso_original[molecule] = {};

        iso_ao[molecule] = {};
        aniso_ao[molecule] = {};

        // Loop through field strengths
        for (double field_mag : field_values) {

            // Original method
            arma::mat alpha_original = calculate_polarizability_tensor(atoms_file_path.string(), p, q, false, field_mag);
            double iso_original_val = isotropic_polarizability(alpha_original);
            double aniso_original_val = anisotropic_polarizability(alpha_original);

            iso_original[molecule].push_back(iso_original_val);
            aniso_original[molecule].push_back(aniso_original_val);

            // Dipole AO integral method
            arma::mat alpha_ao_mat = calculate_polarizability_tensor(atoms_file_path.string(), p, q, true, field_mag);
            double iso_ao_val = isotropic_polarizability(alpha_ao_mat);
            double aniso_ao_val = anisotropic_polarizability(alpha_ao_mat);

            iso_ao[molecule].push_back(iso_ao_val);
            aniso_ao[molecule].push_back(aniso_ao_val);

            std::cout << "Finished " << molecule << " at field " << field_mag << std::endl;
        }
    }

    // Write out the csvs
    write_field_dependence_csv(output_dir + "/iso_polarizability_original_vs_field.csv", field_values, molecule_names, iso_original);
    write_field_dependence_csv(output_dir + "/aniso_polarizability_original_vs_field.csv", field_values, molecule_names, aniso_original);
    write_field_dependence_csv(output_dir + "/iso_polarizability_ao_vs_field.csv", field_values, molecule_names, iso_ao);
    write_field_dependence_csv(output_dir + "/aniso_polarizability_ao_vs_field.csv", field_values, molecule_names, aniso_ao);

    std::cout << "\nField dependence pipeline complete.\n" << "Results written to: " << output_dir << std::endl;

    return 0;
}