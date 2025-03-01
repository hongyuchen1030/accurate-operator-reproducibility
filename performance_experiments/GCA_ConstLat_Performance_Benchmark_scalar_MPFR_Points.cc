#include <tuple>
#include <iostream>
#include <chrono>
#include "simd_fma.hh"
#include <mp++/mp++.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include "GCAConstLatIntersections.h"
#include "scalar_type_benchmark.hh"
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include "scalar_type_benchmark.hh"


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
    static constexpr size_t dataSize = 10000000;
    Eigen::MatrixXd ptsA(dataSize, 3), ptsB(dataSize, 3);
    Eigen::VectorXd lattitudes(dataSize), X(dataSize), Y(dataSize);

    // Define some test points and constant Z using the provided values
    double pointA_x_significand = 5028374390644146;
    int pointA_x_exponent = -53;
    double pointA_y_significand = -7472957205960756;
    int pointA_y_exponent = -53;
    double pointA_z_significand = 6254432438282003;
    int pointA_z_exponent = -79;

    double pointB_x_significand = 5167685454902838;
    int pointB_x_exponent = -53;
    double pointB_y_significand = -7377307466399399;
    int pointB_y_exponent = -53;
    double pointB_z_significand = 4525606513452550;
    int pointB_z_exponent = -78;

    double constZ_significand = 7998403280412384;
    int constZ_exponent = -79;

    // Constructing PointA, PointB, and constZ
    Eigen::Vector3d pointA(
        pointA_x_significand * pow(2.0, pointA_x_exponent),
        pointA_y_significand * pow(2.0, pointA_y_exponent),
        pointA_z_significand * pow(2.0, pointA_z_exponent)
    );

    Eigen::Vector3d pointB(
        pointB_x_significand * pow(2.0, pointB_x_exponent),
        pointB_y_significand * pow(2.0, pointB_y_exponent),
        pointB_z_significand * pow(2.0, pointB_z_exponent)
    );

    double constZ = constZ_significand * pow(2.0, constZ_exponent);

    for (size_t i = 0; i < dataSize; ++i) {
        ptsA.row(i) = pointA;
        ptsB.row(i) = pointB;
    }

    lattitudes.setConstant(constZ);

    // Add small random perturbations to the points and latitudes
    ptsA += 1e-12 * Eigen::MatrixXd::Random(ptsA.rows(), ptsA.cols());
    ptsB += 1e-12 * Eigen::MatrixXd::Random(ptsB.rows(), ptsB.cols());
    lattitudes += 1e-12 * Eigen::VectorXd::Random(lattitudes.size());

    // Reference results for comparison
    Eigen::VectorXd X_true(dataSize), Y_true(dataSize);
    compute_values_EFT<double>(ptsA, ptsB, lattitudes, X_true, Y_true);  // Ensure this function is declared properly

    Eigen::VectorXd X_true_fp(dataSize), Y_true_fp(dataSize);
    compute_values_new<double>(ptsA, ptsB, lattitudes, X_true_fp, Y_true_fp);
    compute_values_old<double>(ptsA, ptsB, lattitudes, X_true_fp, Y_true_fp);

    // Initialize CSV file for output
    std::string output_directory = "../output";  // relative path to the output directory
    std::string output_file = output_directory + "/average_time_scalar_MPFR_points.csv";

    // Ensure the output directory exists
    create_output_directory(output_directory);

    // Open the CSV file for writing
    std::ofstream csv_file(output_file);
    csv_file << std::scientific << std::setprecision(16);

    // Write CSV header
    csv_file << "float_OldEqn,float_NewEqn,float_BaselineEqn,our_method";
    for (int precision = 16; precision <= 32; precision += 16) {
        csv_file << ",MPFR" << precision << "_OldEqn,MPFR" << precision << "_NewEqn,MPFR" << precision << "_BaselineEqn";
    }
    csv_file << ",Quadruple_OldEqn,Quadruple_NewEqn,Quadruple_BaselineEqn\n";

    size_t num_tests = 0, passed = 0;
    
    // Function to calculate and store errors
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

    // Run the benchmark and save the average times
    for (size_t repeat = 0; repeat < 1; ++repeat) {
        std::cout << "\n========= Naive Floating Point Results =========\n" << std::endl;
        auto [time_new, time_old, time_baseline] = run_benchmark<double>(ptsA, ptsB, lattitudes, X, Y);

        double avg_time_new = time_new / (numTests * dataSize);
        double avg_time_old = time_old / (numTests * dataSize);
        double avg_time_baseline = time_baseline / (numTests * dataSize);

        // Print out the FP results
        std::cout<< " \n Naive FP baseline equation: " << avg_time_baseline <<" seconds \n" << std::endl;
        std::cout<< " \n Naive FP old equation: " << avg_time_old <<" seconds \n" << std::endl;
        std::cout<< " \n Naive FP new equation: " << avg_time_new <<" seconds \n" << std::endl;


        // **EFT Results (our_method)**
        std::cout << "\n========= EFT Results =========\n" << std::endl;
        auto eft_duration = run_benchmark_EFT<double>(ptsA, ptsB, lattitudes, X, Y);
        double avg_eft_time = eft_duration / (numTests * dataSize);
        // verify(1);
        X.setZero();
        Y.setZero();
        std::cout<< " \n Our EFT: " << avg_eft_time <<" seconds \n" << std::endl;

        // Write results to CSV
        csv_file << avg_time_old << "," << avg_time_new << "," << avg_time_baseline << "," << avg_eft_time;

        // MPFR test
        for (int precision = 16; precision <= 32; precision += 16) {
            // Print the current precision
            std::cout << "\n========= MPFR precision: " << precision << " bits =========\n" << std::endl;

            auto [time_mpfr_new, time_mpfr_old, time_mpfr_baseline] = handle_mpfr_precision(precision, ptsA, ptsB, lattitudes, X, Y);

            // Write the results to the CSV file
            csv_file
                    << "," << time_mpfr_old / (numTests * dataSize)
                    << "," << time_mpfr_new / (numTests * dataSize)
                    << "," << time_mpfr_baseline / (numTests * dataSize) ;
        }


        std::cout << "\n========= Quadruple precision Results =========\n" << std::endl;

        auto [time_new_quad, time_old_quad, time_baseline_quad] = run_quadruple_benchmark(ptsA, ptsB, lattitudes, X, Y);

        // Write Quadruple results to CSV
        csv_file << "," << time_old_quad / (numTests * dataSize)
                 << "," << time_new_quad / (numTests * dataSize)
                 << "," << time_baseline_quad / (numTests * dataSize)
                 << "\n";
    }

    std::cout << passed << "/" << num_tests  << " tests passed" << std::endl;

    // Close the CSV file
    csv_file.close();

    return 0;
}