#include <tuple>
#include <iostream>
#include <chrono>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>

#include "forward_declarations.hh"

#include <tbb/parallel_for.h>
#include <tbb/global_control.h>
std::unique_ptr<tbb::global_control> g_global_control;
void set_max_num_tbb_threads(int num_threads) {
    if (num_threads < 1) throw std::runtime_error("num_threads must be >= 1");
    g_global_control = std::make_unique<tbb::global_control>(tbb::global_control::parameter::max_allowed_parallelism, num_threads);
}


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
    Eigen::MatrixXd ptsA(dataSize, 3), ptsB(dataSize, 3);
    Eigen::VectorXd lattitudes(dataSize), X(dataSize), Y(dataSize);

    // Define test points and constant Z using the provided values
    Eigen::Vector3d pointA(5028374390644146 * pow(2.0, -53), -7472957205960756 * pow(2.0, -53), 6254432438282003 * pow(2.0, -79));
    Eigen::Vector3d pointB(5167685454902838 * pow(2.0, -53), -7377307466399399 * pow(2.0, -53), 4525606513452550 * pow(2.0, -78));
    double constZ = 7998403280412384 * pow(2.0, -79);

    for (size_t i = 0; i < dataSize; ++i) {
        ptsA.row(i) = pointA;
        ptsB.row(i) = pointB;
    }

    lattitudes.setConstant(constZ);
    ptsA += 1e-12 * Eigen::MatrixXd::Random(ptsA.rows(), ptsA.cols());
    ptsB += 1e-12 * Eigen::MatrixXd::Random(ptsA.rows(), ptsA.cols());
    lattitudes += 1e-12 * Eigen::VectorXd::Random(ptsA.rows());

    Eigen::VectorXd X_true(dataSize), Y_true(dataSize);
    compute_values_EFT<double, true>(ptsA, ptsB, lattitudes, X_true, Y_true);

    Eigen::VectorXd X_true_fp(dataSize), Y_true_fp(dataSize);
    compute_values<double, true>(ptsA, ptsB, lattitudes, X_true_fp, Y_true_fp);

    size_t num_tests = 0, passed = 0;
    auto verify = [&](size_t vector_width) {
        double x_error = (X - X_true).norm();
        double y_error = (Y - Y_true).norm();
        if ((x_error > 0) || (y_error > 0)) {
            std::cerr << "Error with vector width " << vector_width << ": " << x_error << " " << y_error << std::endl;
        } else {
            ++passed;
        }
        ++num_tests;
    };

    auto verify_fp = [&](size_t vector_width) {
        double x_error = (X - X_true_fp).norm();
        double y_error = (Y - Y_true_fp).norm();
        if ((x_error > 0) || (y_error > 0)) {
            std::cerr << "Error with vector width " << vector_width << ": " << x_error << " " << y_error << std::endl;
        } else {
            ++passed;
        }
        ++num_tests;
    };

        // Initialize CSV file for output
    std::string output_directory = "../output";  // relative path to the output directory
    std::string output_file = output_directory + "/average_time_parallel.csv";

    // Ensure the output directory exists
    create_output_directory(output_directory);

    // Open the CSV file for writing
    std::ofstream csv_file(output_file);
    csv_file << std::scientific << std::setprecision(16);

    // Write CSV header
    csv_file << "Threads,Scalar,EFT_SIMD_2,EFT_SIMD_4,EFT_SIMD_8,FP_Scalar,FP_SIMD_2,FP_SIMD_4,FP_SIMD_8\n";

    for (size_t repeat = 0; repeat < 2; ++repeat) {
        // Parallel versions
        for (size_t numThreads = 1; numThreads < 32; numThreads *= 2) {
            std::cout << std::endl;
            set_max_num_tbb_threads(numThreads);

            // Run benchmarks for different vector widths
            double scalar_time_our = run_benchmark_EFT<double, true>(ptsA, ptsB, lattitudes, X, Y) / (numTest_hardcoded * dataSize);
            std::cout << "Parallel (" << numThreads << " threads) Scalar time: " << scalar_time_our << " seconds" << std::endl;
            verify(1); X.setZero(); Y.setZero();

            double vec_width2_time_our = run_benchmark_EFT<Vec<2>, true>(ptsA, ptsB, lattitudes, X, Y) / (numTest_hardcoded * dataSize);
            std::cout << "Parallel (" << numThreads << " threads) SIMD time (width 2): " << vec_width2_time_our << " seconds" << std::endl;
            verify(2); X.setZero(); Y.setZero();

            double vec_width4_time_our = run_benchmark_EFT<Vec<4>, true>(ptsA, ptsB, lattitudes, X, Y) / (numTest_hardcoded * dataSize);
            std::cout << "Parallel (" << numThreads << " threads) SIMD time (width 4): " << vec_width4_time_our << " seconds" << std::endl;
            verify(4); X.setZero(); Y.setZero();

            double vec_width8_time_our = run_benchmark_EFT<Vec<8>, true>(ptsA, ptsB, lattitudes, X, Y) / (numTest_hardcoded * dataSize);
            std::cout << "Parallel (" << numThreads << " threads) SIMD time (width 8): " << vec_width8_time_our << " seconds" << std::endl;
            verify(8); X.setZero(); Y.setZero();

            // Run benchmarks for floating point versions
            double scalar_time_fp = run_benchmark<double, true>(ptsA, ptsB, lattitudes, X, Y) / (numTest_hardcoded * dataSize);
            std::cout << "Parallel (" << numThreads << " threads) FP Scalar time: " << scalar_time_fp << " seconds" << std::endl;
            verify_fp(1); X.setZero(); Y.setZero();

            double vec_width2_time_fp = run_benchmark<Vec<2>, true>(ptsA, ptsB, lattitudes, X, Y) / (numTest_hardcoded * dataSize);
            std::cout << "Parallel (" << numThreads << " threads) FP SIMD time (width 2): " << vec_width2_time_fp << " seconds" << std::endl;
            verify_fp(2); X.setZero(); Y.setZero();

            double vec_width4_time_fp = run_benchmark<Vec<4>, true>(ptsA, ptsB, lattitudes, X, Y) / (numTest_hardcoded * dataSize);
            std::cout << "Parallel (" << numThreads << " threads) FP SIMD time (width 4): " << vec_width4_time_fp << " seconds" << std::endl;
            verify_fp(4); X.setZero(); Y.setZero();

            double vec_width8_time_fp = run_benchmark<Vec<8>, true>(ptsA, ptsB, lattitudes, X, Y) / (numTest_hardcoded * dataSize);
            std::cout << "Parallel (" << numThreads << " threads) FP SIMD time (width 8): " << vec_width8_time_fp << " seconds" << std::endl;
            verify_fp(8); X.setZero(); Y.setZero();

            // Write the results to the CSV file
            csv_file << numThreads << ","
                     << scalar_time_our << ","
                     << vec_width2_time_our << ","
                     << vec_width4_time_our << ","
                     << vec_width8_time_our << ","
                     << scalar_time_fp << ","
                     << vec_width2_time_fp << ","
                     << vec_width4_time_fp << ","
                     << vec_width8_time_fp << "\n";
        }
    }

    // Close the CSV file
    csv_file.close();




    std::cout << passed << "/" << num_tests << " tests passed" << std::endl;
    return 0;
}
