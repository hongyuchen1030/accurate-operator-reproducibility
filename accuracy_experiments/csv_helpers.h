//
// Created by hongyu chen on 4/25/24.
//

#ifndef CODES_CSV_HELPER_H
#define CODES_CSV_HELPER_H
#include <fstream>
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <Eigen/Dense>
#include <boost/math/tools/roots.hpp>
#include <filesystem> // Ensure this is included
#include <mp++/mp++.hpp>
namespace fs = std::filesystem;


using namespace boost::multiprecision;
template<typename T>
using V3_T = Eigen::Matrix<T, 3, 1>;
template<typename T>
using Arc_T = Eigen::Matrix<T, 3, 2>; // Matrix of 3 rows and 2 columns


template<typename T>
struct PrecisionTraits;

template<>
struct PrecisionTraits<float> { static constexpr size_t p = 24; };  // 24 bits for the mantissa

template<>
struct PrecisionTraits<double> { static constexpr size_t p = 53; }; // 53 bits for the mantissa
// Function to log errors
void logError(const std::string& methodName, const std::string& message, double value, int64_t significand, int64_t exponent) {
    std::cerr << "Error in method: " << methodName << ", Value: " << value << ", Message: " << message 
              << ", Significand: " << significand << ", Exponent: " << exponent << std::endl;
}
struct DecomposedFloat {
    long significand;
    long exponent;

    DecomposedFloat() : significand(0), exponent(0) {}

    // Constructor to directly set significand and exponent
    DecomposedFloat(long sig, long exp) : significand(sig), exponent(exp) {}

    template<typename T>
    DecomposedFloat(T v) {
        if (v == 0) {
            significand = 0;
            exponent = 0;
        } else {
            exponent = std::ilogb(v) - (std::numeric_limits<T>::digits - 1);
            T value_exp = std::ldexp(v, -exponent);
            significand = long(value_exp);
            if (significand == 0) {
                exponent = 0;
            }
            if (T(significand) != value_exp) {
                logError("DecomposedFloat", "Integer part was not integer!", v, significand, exponent);
                throw std::runtime_error("Integer part was not integer!");
            }
            if (toFloat<T>() != v) {
                logError("DecomposedFloat", "Could not reconstruct the original value!", v, significand, exponent);
                throw std::runtime_error("Could not reconstruct the original value!");
            }
        }
    }

    // Output stream operator
    friend std::ostream& operator<<(std::ostream& os, const DecomposedFloat &df) {
        return os << df.significand << " * 2**" << df.exponent;
    }

    template<typename T>
    T toFloat() const {
        return std::ldexp(T(significand), exponent);
    }

    // Convenience function to directly output as CSV
    std::string toCSV() const {
        return std::to_string(significand) + "," + std::to_string(exponent);
    }
};

std::string formatLatitude(double latitude) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(6) << latitude;
    std::string latStr = stream.str();
    for (auto& c : latStr) {
        if (c == '.') {
            c = '_';
        }
    }
    return latStr;
}

std::string formatOffset(double offset) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(16) << offset;
    std::string latStr = stream.str();
    for (auto& c : latStr) {
        if (c == '.') {
            c = '_';
        }
    }
    return latStr;
}

std::vector<double> parseDoubles(const std::string& input) {
    std::vector<double> doubles;
    std::stringstream ss(input);
    std::string item;
    while (getline(ss, item, ',')) {
        doubles.push_back(std::atof(item.c_str()));
    }
    return doubles;
}

std::vector<int> parseInts(const std::string& input) {
    std::vector<int> ints;
    std::stringstream ss(input);
    std::string item;
    while (getline(ss, item, ',')) {
        ints.push_back(std::stoi(item));
    }
    return ints;
}



// Function to read great circle arcs from a CSV file into a vector
std::vector<Arc_T<double>> readGreatCircleArcsFromCSV(const std::string& filename) {
    std::vector<Arc_T<double>> arcs;
    std::ifstream csvFile(filename);
    std::string line;
    std::getline(csvFile, line); // Skip the header line

    while (std::getline(csvFile, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<double> values;

        // Extract each comma-separated value
        while (std::getline(ss, value, ',')) {
            values.push_back(std::stod(value));
        }

        // Create an Arc_T and push it to the vector
        Arc_T<double> arc;
        arc << V3_T<double>(values[0], values[1], values[2]),
                V3_T<double>(values[3], values[4], values[5]);
        arcs.push_back(arc);
    }

    csvFile.close();
    return arcs;
}

// Function to generate and write random great circle arcs to a CSV file
static void writeGreatCircleArcsToCSV(
        const std::string& filename,
        std::vector<Arc_T<double>> arcs, const int print_precision = 32) {
    // Check if the arcs vector is empty;
    if (arcs.empty()) {
        std::cerr << "No arcs to write to CSV." << std::endl;
        return;
    }

    // Open file stream to write
    std::ofstream csvFile(filename);

    // Write CSV headers
    csvFile << "\"gca_pointA_x\",\"gca_pointA_y\",\"gca_pointA_z\",\"gca_pointB_x\",\"gca_pointB_y\",\"gca_pointB_z\",\"const_z\"\n";


    // Set the stream to scientific notation
    csvFile << std::scientific;

    // Set the precision using the member function
    csvFile.precision(print_precision);

    // Write arc data to CSV
    for (const auto& arc : arcs) {
        double const_z = 0.5 * (arc.col(0).z() + arc.col(1).z());
        csvFile << arc.col(0).x() << "," << arc.col(0).y() << "," << arc.col(0).z() << ","
                << arc.col(1).x() << "," << arc.col(1).y() << "," << arc.col(1).z() << ","
                << const_z << "\n";
    }

    // Close the file
    csvFile.close();
}

void writeIVC_Double(const std::vector<intermediate_values_compare<double>>& ivcs, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    // Write headers reflecting significand and exponent for each value
    outFile << "nx_significand,nx_exponent,"
            << "nNorm_significand,nNorm_exponent,"
            << "sSquare_significand,sSquare_exponent,"
            << "s_significand,s_exponent,"
            << "znxnz_significand,znxnz_exponent,"
            << "sny_significand,sny_exponent,"
            << "znxnzsny_significand,znxnzsny_exponent,"
            << "nxSquarenySquare_significand,nxSquarenySquare_exponent,"
            << "px_significand,px_exponent\n";

    // Decompose each value and write to file
    for (auto ivc : ivcs) {
        outFile << DecomposedFloat(ivc.nx).toCSV() << ","
                << DecomposedFloat(ivc.nNorm).toCSV() << ","
                << DecomposedFloat(ivc.sSquare).toCSV() << ","
                << DecomposedFloat(ivc.s).toCSV() << ","
                << DecomposedFloat(ivc.znxnz).toCSV() << ","
                << DecomposedFloat(ivc.sny).toCSV() << ","
                << DecomposedFloat(ivc.znxnz_sny).toCSV() << ","
                << DecomposedFloat(ivc.nxSquarenySquare).toCSV() << ","
                << DecomposedFloat(ivc.px).toCSV() << "\n";

    }


    outFile.close();
}



void writeIVC_Quadruple(const std::vector<intermediate_values_compare<mppp::real128>>& ivcs_quadruple, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    // Write headers with precision info
    outFile << "nx,"
            << "nNorm,"
            << "sSquare,"
            << "s,"
            << "znxnz,"
            << "sny,"
            << "znxnzsny,"
            << "nxSquarenySquare,"
            << "px\n";

    // Set stream to scientific notation once
    outFile.precision(34); // 34 is typically the precision for quadruple precision (128-bit)
    outFile << std::scientific;

    for (const auto& iv : ivcs_quadruple) {
        outFile << iv.nx << ","
                << iv.nNorm << ","
                << iv.sSquare << ","
                << iv.s << ","
                << iv.znxnz << ","
                << iv.sny << ","
                << iv.znxnz_sny << ","
                << iv.nxSquarenySquare << ","
                << iv.px << "\n";
    }

    outFile.close();
}


void writeIVC_MPFR(const std::vector<intermediate_values_compare<mpfr_float>>& ivcs, const std::string& filename, int precision) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }


    // Write headers with precision info
    outFile << "nx,"
            << "nNorm,"
            << "sSquare,"
            << "s,"
            << "znxnz,"
            << "sny,"
            << "znxnzsny,"
            << "nxSquarenySquare,"
            << "px\n";



    // Set stream to scientific notation once
    outFile.precision(precision);

    for (auto iv : ivcs) {
//        std::cout << "MPFR precision: " << precision<<" Has Px: " <<iv.nx <<std::endl;
        outFile<< std::scientific << iv.nx << ","
                << iv.nNorm << ","
                << iv.sSquare << ","
                << iv.s << ","
                << iv.znxnz << ","
                << iv.sny << ","
                << iv.znxnz_sny << ","
                << iv.nxSquarenySquare << ","
                << iv.px << "\n";

    }

    outFile.close();
}


void writeIV_Our(const std::vector<intermediate_values_our>& iv_our, const std::string& filename, const int print_precision = 32) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    // Write headers for tuples
    outFile << "nx_1_significand,nx_1_exponent,nx_2_significand,nx_2_exponent,"
            << "nNorm_1_significand,nNorm_1_exponent,nNorm_2_significand,nNorm_2_exponent,"
            << "sSquare_1_significand,sSquare_1_exponent,sSquare_2_significand,sSquare_2_exponent,"
            << "s_1_significand,s_1_exponent,s_2_significand,s_2_exponent,"
            << "znxnz_1_significand,znxnz_1_exponent,znxnz_2_significand,znxnz_2_exponent,"
            << "sny_1_significand,sny_1_exponent,sny_2_significand,sny_2_exponent,"
            << "znxnzsny_significand,znxnzsny_exponent,"
            << "nxSquarenySquare_significand,nxSquarenySquare_exponent,"
            << "px_significand,px_exponent\n";


    // Write each entry in iv_our
    for (const auto& iv : iv_our) {
        outFile << DecomposedFloat(std::get<0>(iv.nx_pair)).toCSV() << ","
                << DecomposedFloat(std::get<1>(iv.nx_pair)).toCSV() << ","
                << DecomposedFloat(std::get<0>(iv.nNorm_pair)).toCSV() << ","
                << DecomposedFloat(std::get<1>(iv.nNorm_pair)).toCSV() << ","
                << DecomposedFloat(std::get<0>(iv.sSquare_pair)).toCSV() << ","
                << DecomposedFloat(std::get<1>(iv.sSquare_pair)).toCSV() << ","
                << DecomposedFloat(std::get<0>(iv.s_pair)).toCSV() << ","
                << DecomposedFloat(std::get<1>(iv.s_pair)).toCSV() << ","
                << DecomposedFloat(std::get<0>(iv.znxnz_pair)).toCSV() << ","
                << DecomposedFloat(std::get<1>(iv.znxnz_pair)).toCSV() << ","
                << DecomposedFloat(std::get<0>(iv.sny_pair)).toCSV() << ","
                << DecomposedFloat(std::get<1>(iv.sny_pair)).toCSV() << ","
                << DecomposedFloat(iv.znxnz_sny).toCSV() << ","
                << DecomposedFloat(iv.nxSquarenySquare).toCSV() << ","
                << DecomposedFloat(iv.px).toCSV() << "\n";
    }

    outFile.close();

}


// Function to generate and write random great circle arcs to a CSV file in an exponential format
static void WriteGreatCircleArcsToCSVExponent(
        const std::string& filename,
        std::vector<Arc_T<double>> arcs, std::vector<double>& constZs) {
    // Check if the arcs vector is empty;
    if (arcs.empty()) {
        std::cerr << "No arcs to write to CSV." << std::endl;
        return;
    }

    // Check if the file exists and delete it
    if (fs::exists(filename)) {
        try {
            fs::remove(filename);
        } catch (const fs::filesystem_error& e) {
            std::cerr << "Error deleting file: " << e.what() << std::endl;
            return; // Return with error code
        }
    }

    // Open file stream to write
    std::ofstream csvFile(filename);

    // Write CSV headers
    csvFile << "pointA_x_significand,pointA_x_exponent,pointA_y_significand,pointA_y_exponent,pointA_z_significand,pointA_z_exponent,"
               "pointB_x_significand,pointB_x_exponent,pointB_y_significand,pointB_y_exponent,pointB_z_significand,pointB_z_exponent,"
               "constZ_significand,constZ_exponent\n";

    // Process each arc and write decomposed values
    for (int i = 0; i < arcs.size(); ++i) {
        auto arc = arcs[i];
        // Calculate constZ as average of z coordinates of point A and B
        double constZ = constZs[i];  // z coordinates of point A and B
        DecomposedFloat df_constZ(constZ);

        for (int col = 0; col < 2; ++col) {  // Columns for Point A and Point B
            for (int row = 0; row < 3; ++row) {  // Rows for x, y, z
                DecomposedFloat df(arc(row, col));
                csvFile << df.significand << ',' << df.exponent << ',';
            }
        }

        // Write constZ decomposed to file
        csvFile << df_constZ.significand << ',' << df_constZ.exponent << '\n';
    }

    // Close the file
    csvFile.close();
}


// Helper function to reconstruct double from significand and exponent
double reconstructDouble(long significand, long exponent) {
    return std::ldexp(static_cast<double>(significand), exponent);
}

// Function to read great circle arcs from a CSV file into a vector
std::vector<Arc_T<double>> ReadGreatCircleArcsFromCSVExponent(const std::string& filename, std::vector<double>& constZs) {
    std::vector<Arc_T<double>> arcs;
    std::ifstream csvFile(filename);
    if (!csvFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return arcs;
    }

    std::string line;
    std::getline(csvFile, line); // Skip the header line

    while (std::getline(csvFile, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<long> significands;
        std::vector<int> exponents;

        // Extract each comma-separated value and separate significands and exponents
        while (std::getline(ss, value, ',')) {
            long significand;
            int exponent;
            static bool toggle = true; // Toggle between reading significand and exponent

            if (toggle) { // Read significand
                significand = std::stol(value);
                significands.push_back(significand);
            } else { // Read exponent
                exponent = std::stoi(value);
                exponents.push_back(exponent);
            }
            toggle = !toggle; // Flip toggle for next read
        }


        // Reconstruct doubles and create an Arc_T
        Arc_T<double> arc;
        for (int i = 0; i < 6; ++i) {
            int row = i % 3;
            int col = i / 3;
            arc(row, col)= reconstructDouble(significands[i], exponents[i]);
        }
        arcs.push_back(arc);

        // Reconstruct constZ and add to the vector
        double constZ = reconstructDouble(significands[6], exponents[6]);
        constZs.push_back(constZ);
    }


    csvFile.close();
    return arcs;
}

// Function to read and reconstruct data from a CSV file
std::vector<std::vector<double>> ReadAndReconstructArcData(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<double>> data;
    std::string line, cell;

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return data;
    }

    // Skip the header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::vector<double> reconstructedRow;
        int index = 0;
        long significand;
        int exponent;

        while (std::getline(lineStream, cell, ',')) {
            if (index % 2 == 0) {  // Significand
                significand = std::stol(cell);
            } else {              // Exponent
                exponent = std::stoi(cell);
                // Reconstruct the double and add to the row
                reconstructedRow.push_back(reconstructDouble(significand, exponent));
            }
            index++;
        }
        data.push_back(reconstructedRow);
    }

    file.close();
    return data;
}

// Function to read directly real number data from the reconstructed CSV file
std::vector<std::vector<double>> ReadDirectDoubleData(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::vector<double>> data;
    std::string line, cell;

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return data;
    }

    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::vector<double> row;
        while (std::getline(lineStream, cell, ',')) {
            row.push_back(std::stod(cell));
        }
        data.push_back(row);
    }

    file.close();
    return data;
}



#endif //CODES_CSV_HELPER_H
