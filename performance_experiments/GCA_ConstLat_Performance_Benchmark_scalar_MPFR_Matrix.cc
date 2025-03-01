#include <tuple>
#include <iostream>
#include <chrono>
#include "simd_fma.hh"
#include <mp++/mp++.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include "GCAConstLatIntersections.h"
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>


template<typename T, size_t N>
using VN_T = Eigen::Matrix<T, N, 1>;
template<typename T>
using V3_T = VN_T<T, 3>;


using namespace boost::multiprecision;
template<unsigned Precision>
using mpfr_static_t = number<mpfr_float_backend<Precision>, et_on>;



template <typename T>
constexpr size_t vecWidth() {
    if constexpr (std::is_arithmetic_v<T>) return 1;
    else return T::RowsAtCompileTime;
}

template<size_t VEC_WIDTH>
using Vec = Eigen::Array<double, VEC_WIDTH, 1>;

constexpr size_t numTests = 100;

template <typename T>
void compute_values_EFT(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    const size_t size = ptsA.rows();
    constexpr size_t stride = 1;


    std::array<T, 3> x1, x2;
    T z0;
    for (size_t i = 0; i < size; i += stride) {

        for (size_t c = 0; c < 3; ++c) {
            x1[c] = ptsA(i, c);
            x2[c] = ptsB(i, c);
        }
        z0 = lattitudes[i];


        auto [px, py] = AccurateIntersection::gca_constLat_intersection_accurate_performance(x1, x2, z0);

        X[i] = px;
        Y[i] = py;

    }
}

template <typename T>
double run_benchmark_EFT(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < numTests; ++i) {
         compute_values_EFT<T>(ptsA, ptsB, lattitudes, X, Y);
    }

    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;

    return elapsed.count();
}


template <typename T>
void compute_values_new(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &ptsA, 
                    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &ptsB, 
                    const Eigen::Matrix<T, Eigen::Dynamic, 1> &lattitudes, 
                    Eigen::Matrix<T, Eigen::Dynamic, 1> &X, 
                    Eigen::Matrix<T, Eigen::Dynamic, 1> &Y) {
    const size_t size = ptsA.rows();
    constexpr size_t stride = 1;

    V3_T<T> x1, x2;
    T z0;
    
    for (size_t i = 0; i < size; i += stride) {
        for (size_t c = 0; c < 3; ++c) {
            x1(c) = ptsA(i, c);  // Assign directly for scalar T
            x2(c) = ptsB(i, c);
        }
        z0 = lattitudes(i);

        // Call the provided Intersection Function
        auto [px, py] = gca_constLat_intersection_coordinates_newEqn_performance(x1, x2, z0);

        X[i] = px;
        Y[i] = py;
    }
}

template <typename T>
void compute_values_old(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &ptsA, 
                    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &ptsB, 
                    const Eigen::Matrix<T, Eigen::Dynamic, 1> &lattitudes, 
                    Eigen::Matrix<T, Eigen::Dynamic, 1> &X, 
                    Eigen::Matrix<T, Eigen::Dynamic, 1> &Y) {
    const size_t size = ptsA.rows();
    constexpr size_t stride = 1;

    V3_T<T> x1, x2;
    T z0;
    
    for (size_t i = 0; i < size; i += stride) {
        for (size_t c = 0; c < 3; ++c) {
            x1(c) = ptsA(i, c);  // Assign directly for scalar T
            x2(c) = ptsB(i, c);
        }
        z0 = lattitudes(i);

        // Call the provided Intersection Function
        auto [px, py] = gca_constLat_intersection_coordinates_oldEqn_performance(x1, x2, z0);

        X[i] = px;
        Y[i] = py;
    }
}

template <typename T>
void compute_values_baseline(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &ptsA, 
                    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &ptsB, 
                    const Eigen::Matrix<T, Eigen::Dynamic, 1> &lattitudes, 
                    Eigen::Matrix<T, Eigen::Dynamic, 1> &X, 
                    Eigen::Matrix<T, Eigen::Dynamic, 1> &Y) {
    const size_t size = ptsA.rows();
    constexpr size_t stride = 1;

    V3_T<T> x1, x2;
    T z0;
    
    for (size_t i = 0; i < size; i += stride) {
        for (size_t c = 0; c < 3; ++c) {
            x1(c) = ptsA(i, c);  // Assign directly for scalar T
            x2(c) = ptsB(i, c);
        }
        z0 = lattitudes(i);

        // Call the provided Intersection Function
        auto [px, py] = gca_constLat_intersection_coordinates_baselineEqn_performance(x1, x2, z0);

        X[i] = px;
        Y[i] = py;
    }
}



template <typename T>
std::tuple<double, double, double> run_benchmark(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &ptsA, 
                                                 const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &ptsB, 
                                                 const Eigen::Matrix<T, Eigen::Dynamic, 1> &lattitudes, 
                                                 Eigen::Matrix<T, Eigen::Dynamic, 1> &X, 
                                                 Eigen::Matrix<T, Eigen::Dynamic, 1> &Y) {
    // Get the intersection functions for the three different equations (new, old, baseline)

    // To store the timings for each method
    double time_new = 0.0, time_old = 0.0, time_baseline = 0.0;

    // Measure time for the "new" method
    {
        auto start_new = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < numTests; ++i) {
            // Do not use std::move on X and Y here since we need to reuse them in the loop
            compute_values_new<T>(ptsA, ptsB, lattitudes, X, Y);
        }
        std::chrono::duration<double> elapsed_new = std::chrono::high_resolution_clock::now() - start_new;
        time_new = elapsed_new.count();
    }

    // Measure time for the "old" method
    {
        auto start_old = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < numTests; ++i) {
            compute_values_old<T>(ptsA, ptsB, lattitudes, X, Y);
        }
        std::chrono::duration<double> elapsed_old = std::chrono::high_resolution_clock::now() - start_old;
        time_old = elapsed_old.count();
    }

    // Measure time for the "baseline" method
    {
        auto start_baseline = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < numTests; ++i) {
            compute_values_baseline<T>(ptsA, ptsB, lattitudes, X, Y);
        }
        std::chrono::duration<double> elapsed_baseline = std::chrono::high_resolution_clock::now() - start_baseline;
        time_baseline = elapsed_baseline.count();
    }

    // Return the results as a tuple
    return {time_new, time_old, time_baseline};
}



// Helper function to run MPFR benchmark based on precision
// Helper function to run MPFR benchmark based on precision
template<size_t Precision>
std::tuple<double, double, double> run_mpfr_benchmark(const Eigen::MatrixXd &ptsA, 
                                                      const Eigen::MatrixXd &ptsB, 
                                                      const Eigen::VectorXd &lattitudes) {
    size_t dataSize = ptsA.rows(); // Use runtime data size
    std::cout << "\n========= MPFR Results for Precision " << Precision << " =========\n" << std::endl;

    // Convert ptsA, ptsB, and lattitudes to MPFR format with dynamic size
    Eigen::Matrix<mpfr_static_t<Precision>, Eigen::Dynamic, 3> ptsA_mpfr(dataSize, 3), ptsB_mpfr(dataSize, 3);
    Eigen::Matrix<mpfr_static_t<Precision>, Eigen::Dynamic, 1> latitudes_mpfr(dataSize, 1);

    for (size_t i = 0; i < dataSize; ++i) {
        ptsA_mpfr.row(i) = ptsA.row(i).template cast<mpfr_static_t<Precision>>();
        ptsB_mpfr.row(i) = ptsB.row(i).template cast<mpfr_static_t<Precision>>();
        latitudes_mpfr(i) = mpfr_static_t<Precision>(lattitudes(i));
    }

    // Prepare MPFR X and Y matrices with dynamic size
    Eigen::Matrix<mpfr_static_t<Precision>, Eigen::Dynamic, 1> X_mpfr(dataSize, 1), Y_mpfr(dataSize, 1);
    
    // Run the benchmark and get the timings for MPFR
    auto [time_new_mpfr, time_old_mpfr, time_baseline_mpfr] = run_benchmark<mpfr_static_t<Precision>>(ptsA_mpfr, ptsB_mpfr, latitudes_mpfr, X_mpfr, Y_mpfr);

    // // Print MPFR results
    // std::cout << "New equation MPFR time (Precision " << Precision << "): " << time_new_mpfr << " seconds" << std::endl;
    // std::cout << "Old equation MPFR time (Precision " << Precision << "): " << time_old_mpfr << " seconds" << std::endl;
    // std::cout << "Baseline equation MPFR time (Precision " << Precision << "): " << time_baseline_mpfr << " seconds" << std::endl;

    // Return the timings as a tuple
    return {time_new_mpfr, time_old_mpfr, time_baseline_mpfr};
}

// Function to handle switch-case for MPFR precision without known dataSize
std::tuple<double, double, double> handle_mpfr_precision(int precision, const Eigen::MatrixXd &ptsA, 
                                                         const Eigen::MatrixXd &ptsB, 
                                                         const Eigen::VectorXd &lattitudes) {
    switch (precision) {
        case 16: return run_mpfr_benchmark<16>(ptsA, ptsB, lattitudes);
        case 17: return run_mpfr_benchmark<17>(ptsA, ptsB, lattitudes);
        case 18: return run_mpfr_benchmark<18>(ptsA, ptsB, lattitudes);
        case 19: return run_mpfr_benchmark<19>(ptsA, ptsB, lattitudes);
        case 20: return run_mpfr_benchmark<20>(ptsA, ptsB, lattitudes);
        case 21: return run_mpfr_benchmark<21>(ptsA, ptsB, lattitudes);
        case 22: return run_mpfr_benchmark<22>(ptsA, ptsB, lattitudes);
        case 23: return run_mpfr_benchmark<23>(ptsA, ptsB, lattitudes);
        case 24: return run_mpfr_benchmark<24>(ptsA, ptsB, lattitudes);
        case 25: return run_mpfr_benchmark<25>(ptsA, ptsB, lattitudes);
        case 26: return run_mpfr_benchmark<26>(ptsA, ptsB, lattitudes);
        case 27: return run_mpfr_benchmark<27>(ptsA, ptsB, lattitudes);
        case 28: return run_mpfr_benchmark<28>(ptsA, ptsB, lattitudes);
        case 29: return run_mpfr_benchmark<29>(ptsA, ptsB, lattitudes);
        case 30: return run_mpfr_benchmark<30>(ptsA, ptsB, lattitudes);
        case 31: return run_mpfr_benchmark<31>(ptsA, ptsB, lattitudes);
        case 32: return run_mpfr_benchmark<32>(ptsA, ptsB, lattitudes);
        default:
            std::cerr << "Unsupported MPFR precision: " << precision << std::endl;
            return {0.0, 0.0, 0.0}; // Return a default tuple in case of an error
    }
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
    static constexpr size_t dataSize = 1000000;
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

    // Initialize CSV file for output
    std::string output_directory = "../output";  // relative path to the output directory
    std::string output_file = output_directory + "/average_time.csv";

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
    for (size_t repeat = 0; repeat < 2; ++repeat) {
        std::cout << "\n========= Naive Floating Point Results =========\n" << std::endl;
        auto [time_new, time_old, time_baseline] = run_benchmark<double>(ptsA, ptsB, lattitudes, X, Y);

        double avg_time_new = time_new / (numTests * dataSize);
        double avg_time_old = time_old / (numTests * dataSize);
        double avg_time_baseline = time_baseline / (numTests * dataSize);

        // Print out the FP results
        std::cout<< " \n Naive FP new equation: " << avg_time_new <<" seconds \n" << std::endl;


        // **EFT Results (our_method)**
        std::cout << "\n========= EFT Results =========\n" << std::endl;
        auto eft_duration = run_benchmark_EFT<double>(ptsA, ptsB, lattitudes, X, Y);
        double avg_eft_time = eft_duration / (numTests * dataSize);
        verify(1);
        X.setZero();
        Y.setZero();
        std::cout<< " \n Our EFT: " << avg_eft_time <<" seconds \n" << std::endl;

        // Write results to CSV
        csv_file << avg_time_old << "," << avg_time_new << "," << avg_time_baseline << "," << avg_eft_time;

        // MPFR test
    //     for (int precision = 16; precision <= 32; precision += 16) {
    //         auto [time_mpfr_new, time_mpfr_old, time_mpfr_baseline] = handle_mpfr_precision(precision, ptsA, ptsB, lattitudes);

    //         csv_file << "," << time_mpfr_old / (numTests * dataSize)
    //                  << "," << time_mpfr_new / (numTests * dataSize)
    //                  << "," << time_mpfr_baseline / (numTests * dataSize);
    //     }

    //     std::cout << "\n========= Quadruple precision Results =========\n" << std::endl;

    //     // Quadruple precision test
    //     Eigen::Matrix<mppp::real128, Eigen::Dynamic, 3> ptsA_quad(dataSize, 3), ptsB_quad(dataSize, 3);
    //     Eigen::Matrix<mppp::real128, Eigen::Dynamic, 1> latitudes_quad(dataSize);

    //     for (size_t i = 0; i < dataSize; ++i) {
    //         for (int j = 0; j < 3; ++j) {
    //             ptsA_quad(i, j) = mppp::real128(ptsA(i, j));
    //             ptsB_quad(i, j) = mppp::real128(ptsB(i, j));
    //         }
    //         latitudes_quad(i) = mppp::real128(lattitudes(i));
    //     }

    //     Eigen::Matrix<mppp::real128, Eigen::Dynamic, 1> X_quad(dataSize), Y_quad(dataSize);

    //     auto [time_new_quad, time_old_quad, time_baseline_quad] = run_benchmark<mppp::real128>(ptsA_quad, ptsB_quad, latitudes_quad, X_quad, Y_quad);

    //     // Write Quadruple results to CSV
    //     csv_file << "," << time_old_quad / (numTests * dataSize)
    //              << "," << time_new_quad / (numTests * dataSize)
    //              << "," << time_baseline_quad / (numTests * dataSize)
    //              << "\n";
    }

    std::cout << passed << "/" << num_tests  << " tests passed" << std::endl;

    // Close the CSV file
    csv_file.close();

    return 0;
}