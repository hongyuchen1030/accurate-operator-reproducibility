#include <tuple>
#include <iostream>
#include <chrono>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <Eigen/Dense>
#include <iomanip>

#include "forward_declarations.hh"
#include "generate_input.hh"

// Helper function to create directory if it doesn't exist
void create_output_directory(const std::string& directory) {
    struct stat info;
    if (stat(directory.c_str(), &info) != 0) {
        // Directory does not exist, create it
        if (mkdir(directory.c_str(), 0777) == -1) {
            std::cerr << "Error creating directory: " << directory << std::endl;
        }
    }
}


int main(int argc, const char *argv[]) {
    static constexpr size_t numTest_hardcoded = 100;
    static constexpr size_t dataSize = 100000000;
    // static constexpr size_t dataSize = 16;

    Eigen::MatrixXd ptsA, ptsB;
    Eigen::VectorXd lattitudes;
    generateInput(dataSize, ptsA, ptsB, lattitudes);

    // Storage for computed result.
    Eigen::VectorXd X(dataSize), Y(dataSize);

    Eigen::VectorXd X_true(dataSize), Y_true(dataSize);
    compute_values_EFT<double>(ptsA, ptsB, lattitudes, X_true, Y_true);

    std::cout.precision(18);
    std::cout << "      EFT result: " << X_true[0] << ", " << Y_true[0] << ", " << lattitudes[0] << std::endl;

    Eigen::VectorXd X_true_fp(dataSize), Y_true_fp(dataSize);
    compute_values<double>(ptsA, ptsB, lattitudes, X_true_fp, Y_true_fp);

    std::cout << "Native FP result: " << X_true_fp[0] << ", " << Y_true_fp[0] << ", " << lattitudes[0] << std::endl;
    std::cout << std::endl;

    size_t num_tests = 0, passed = 0;
    auto verify = [&](size_t vector_width) {
        double x_error = (X - X_true).norm();
        double y_error = (Y - Y_true).norm();
        if ((x_error > 0) || (y_error > 0)) {
            std::cerr << "Error with vector width " << vector_width << ": " << x_error << " " << y_error << std::endl;
        }
        else ++passed;
        ++num_tests;
    };


    auto verify_fp = [&](size_t vector_width) {
        double x_error = (X - X_true_fp).norm();
        double y_error = (Y - Y_true_fp).norm();
        if ((x_error > 0) || (y_error > 0)) {
            std::cerr << "Error with vector width " << vector_width << ": " << x_error << " " << y_error << std::endl;
        }
        else ++passed;
        ++num_tests;
    };

    // Initialize CSV file for output
    std::string output_directory = "../output";  // relative path to the output directory
    std::string output_file = output_directory + "/average_time_vectorized.csv";

    // Ensure the output directory exists
    create_output_directory(output_directory);

    // Open the CSV file for writing
    std::ofstream csv_file(output_file);
    csv_file << std::scientific << std::setprecision(16);

    // Write CSV header 
    csv_file << "our_Scalar,our_VecWidth2,our_VecWidth4,our_VecWidth8,FP_Scalar,FP_VecWidth2,FP_VecWidth4,FP_VecWidth8\n";

    for (size_t repeat = 0; repeat < 2; ++repeat) {
        // Run benchmarks and store the results
        double scalar_time_our = run_benchmark_EFT<double, false>(ptsA, ptsB, lattitudes, X, Y);
        verify(1); X.setZero(); Y.setZero();

        double vec_width2_time_our = run_benchmark_EFT<Vec<2>>(ptsA, ptsB, lattitudes, X, Y);
        verify(2); X.setZero(); Y.setZero();

        double vec_width4_time_our = run_benchmark_EFT<Vec<4>>(ptsA, ptsB, lattitudes, X, Y);
        verify(4); X.setZero(); Y.setZero();

        double vec_width8_time_our = run_benchmark_EFT<Vec<8>>(ptsA, ptsB, lattitudes, X, Y);
        verify(8); X.setZero(); Y.setZero();

        double scalar_time_fp = run_benchmark<double>(ptsA, ptsB, lattitudes, X, Y);
        verify_fp(1); X.setZero(); Y.setZero();

        double vec_width2_time_fp = run_benchmark<Vec<2>>(ptsA, ptsB, lattitudes, X, Y);
        verify_fp(2); X.setZero(); Y.setZero();

        double vec_width4_time_fp = run_benchmark<Vec<4>>(ptsA, ptsB, lattitudes, X, Y);
        verify_fp(4); X.setZero(); Y.setZero();

        double vec_width8_time_fp = run_benchmark<Vec<8>>(ptsA, ptsB, lattitudes, X, Y);
        verify_fp(8); X.setZero(); Y.setZero();

        // Print results to the console
        std::cout << "EFT Scalar time: " << scalar_time_our /  (numTest_hardcoded * dataSize)<< " seconds" << std::endl;
        std::cout << "EFT SIMD time (width 2): " << vec_width2_time_our  /  (numTest_hardcoded * dataSize)<< " seconds" << std::endl;
        std::cout << "EFT SIMD time (width 4): " << vec_width4_time_our  /  (numTest_hardcoded * dataSize)<< " seconds" << std::endl;
        std::cout << "EFT SIMD time (width 8): " << vec_width8_time_our  /  (numTest_hardcoded * dataSize)<< " seconds" << std::endl;
        std::cout << "Naive FP Scalar time: " << scalar_time_fp  /  (numTest_hardcoded * dataSize)<< " seconds" << std::endl;
        std::cout << "Naive FP SIMD time (width 2): " << vec_width2_time_fp  /  (numTest_hardcoded * dataSize)<< " seconds" << std::endl;
        std::cout << "Naive FP SIMD time (width 4): " << vec_width4_time_fp /  (numTest_hardcoded * dataSize) << " seconds" << std::endl;
        std::cout << "Naive FP SIMD time (width 8): " << vec_width8_time_fp /  (numTest_hardcoded * dataSize) << " seconds" << std::endl;

        // Write the results to the CSV file
        csv_file << scalar_time_our << "," 
             << vec_width2_time_our << "," 
             << vec_width4_time_our << "," 
             << vec_width8_time_our << "," 
             << scalar_time_fp << "," 
             << vec_width2_time_fp << "," 
             << vec_width4_time_fp << "," 
             << vec_width8_time_fp << "\n";

    }

    // Close the CSV file
    csv_file.close();

    std::cout << passed << "/" << num_tests << " tests passed" << std::endl;

    return 0;
}
