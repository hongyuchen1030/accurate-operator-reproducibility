#ifndef SCALAR_TYPE_BENCHMARK_HH
#define SCALAR_TYPE_BENCHMARK_HH

#include <sys/stat.h>
#include <sys/types.h>
#include <mp++/mp++.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include "GCAConstLatIntersections.h"
#include <Eigen/Dense>
#include <tbb/parallel_for.h>
constexpr size_t numTests = 100;

using namespace boost::multiprecision;
template<unsigned Precision>
using mpfr_static_t = number<mpfr_float_backend<Precision>, et_on>;


template <typename T>
constexpr size_t vecWidth() {
    if constexpr (std::is_arithmetic_v<T>) return 1;
    else return T::RowsAtCompileTime;
}

template <typename T, bool Parallel = false>
void compute_values_EFT(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    const size_t size = ptsA.rows();
    constexpr size_t stride = vecWidth<T>();
    assert((size % stride == 0) && "Size must be a multiple of the vector width");

    if ((size != ptsB.rows()) || (size != lattitudes.rows()) || (size != X.rows()) || (size != Y.rows())) {
        throw std::runtime_error("All passed arrays must have the same size");
    }

    auto processChunk_EFT = [&](size_t chunk_i) {
        const size_t i = chunk_i * stride;
        std::array<T, 3> x1, x2;
        T z0;
        if constexpr (std::is_arithmetic_v<T>) {
            for (size_t c = 0; c < 3; ++c) {
                x1[c] = ptsA(i, c);
                x2[c] = ptsB(i, c);
            }
            z0 = lattitudes[i];
        } else {
            for (size_t c = 0; c < 3; ++c) {
                x1[c] = ptsA.template block<stride, 1>(i, c);
                x2[c] = ptsB.template block<stride, 1>(i, c);
            }
            z0 = lattitudes.template segment<stride>(i);
        }

        auto [px, py] = AccurateIntersection::gca_constLat_intersection_accurate_performance(x1, x2, z0);

        if constexpr (std::is_arithmetic_v<T>) {
            X[i] = px;
            Y[i] = py;
        } else {
            X.template segment<stride>(i) = px;
            Y.template segment<stride>(i) = py;
        }
    };


    if constexpr (Parallel) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size / stride),
                [&](const tbb::blocked_range<size_t> r) {
                    for (size_t i = r.begin(); i != r.end(); ++i)
                        processChunk_EFT(i);
                });
    } else{
        for (size_t i = 0; i < size / stride; ++i) { processChunk_EFT(i); }
    }


}

template <typename T, bool Parallel = false>
double run_benchmark_EFT(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    auto start = std::chrono::high_resolution_clock::now();

    constexpr size_t numTests = 100;
    for (size_t i = 0; i < numTests; ++i) {
         compute_values_EFT<T, Parallel>(ptsA, ptsB, lattitudes, X, Y);
    }

    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;

    return elapsed.count();
}

// New equation performance with casting
template <typename T>
void compute_values_new_cast(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    const size_t size = ptsA.rows();
    V3_T<T> x1, x2;
    T z0;

    for (size_t i = 0; i < size; ++i) {
        for (size_t c = 0; c < 3; ++c) {
            x1(c) = static_cast<T>(ptsA(i, c)); 
            x2(c) = static_cast<T>(ptsB(i, c)); 
        }
        z0 = static_cast<T>(lattitudes(i));
        auto [px, py] = gca_constLat_intersection_coordinates_newEqn_performance(x1, x2, z0);
        X[i] = static_cast<double>(px);
        Y[i] = static_cast<double>(py);
    }
}

// Old equation performance with casting
template <typename T>
void compute_values_old_cast(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    const size_t size = ptsA.rows();
    V3_T<T> x1, x2;
    T z0;

    for (size_t i = 0; i < size; ++i) {
        for (size_t c = 0; c < 3; ++c) {
            x1(c) = static_cast<T>(ptsA(i, c)); 
            x2(c) = static_cast<T>(ptsB(i, c)); 
        }
        z0 = static_cast<T>(lattitudes(i));
        auto [px, py] = gca_constLat_intersection_coordinates_oldEqn_performance(x1, x2, z0);
        X[i] = static_cast<double>(px);
        Y[i] = static_cast<double>(py);
    }
}

// Baseline equation performance with casting
template <typename T>
void compute_values_baseline_cast(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    const size_t size = ptsA.rows();
    V3_T<T> x1, x2;
    T z0;

    for (size_t i = 0; i < size; ++i) {
        for (size_t c = 0; c < 3; ++c) {
            x1(c) = static_cast<T>(ptsA(i, c)); 
            x2(c) = static_cast<T>(ptsB(i, c)); 
        }
        z0 = static_cast<T>(lattitudes(i));
        auto [px, py] = gca_constLat_intersection_coordinates_baselineEqn_performance(x1, x2, z0);
        X[i] = static_cast<double>(px);
        Y[i] = static_cast<double>(py);
    }
}


template <typename T>
void compute_values_new(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
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
void compute_values_old(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
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
void compute_values_baseline(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
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
std::tuple<double, double, double> run_benchmark(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    // Get the intersection functions for the three different equations (new, old, baseline)

    // To store the timings for each method
    double time_new = 0.0, time_old = 0.0, time_baseline = 0.0;

        // Measure time for the "baseline" method
    {
        auto start_baseline = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < numTests; ++i) {
            compute_values_baseline<T>(ptsA, ptsB, lattitudes, X, Y);
        }
        std::chrono::duration<double> elapsed_baseline = std::chrono::high_resolution_clock::now() - start_baseline;
        time_baseline = elapsed_baseline.count();
    }

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



    // Return the results as a tuple
    return {time_new, time_old, time_baseline};
}



// Helper function to run MPFR benchmark based on precision
template<size_t Precision>
std::tuple<double, double, double> run_mpfr_benchmark(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    // Get the intersection functions for the three different equations (new, old, baseline)

    // To store the timings for each method
    double time_new = 0.0, time_old = 0.0, time_baseline = 0.0;

    // Measure time for the "new" method
    {
        auto start_new = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < numTests; ++i) {
            // Do not use std::move on X and Y here since we need to reuse them in the loop
            compute_values_new_cast<mpfr_static_t<Precision>>(ptsA, ptsB, lattitudes, X, Y);
        }
        std::chrono::duration<double> elapsed_new = std::chrono::high_resolution_clock::now() - start_new;
        time_new = elapsed_new.count();
    }

    // Measure time for the "old" method
    {
        auto start_old = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < numTests; ++i) {
            compute_values_old_cast<mpfr_static_t<Precision>>(ptsA, ptsB, lattitudes, X, Y);
        }
        std::chrono::duration<double> elapsed_old = std::chrono::high_resolution_clock::now() - start_old;
        time_old = elapsed_old.count();
    }

    // Measure time for the "baseline" method
    {
        auto start_baseline = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < numTests; ++i) {
            compute_values_baseline_cast<mpfr_static_t<Precision>>(ptsA, ptsB, lattitudes, X, Y);
        }
        std::chrono::duration<double> elapsed_baseline = std::chrono::high_resolution_clock::now() - start_baseline;
        time_baseline = elapsed_baseline.count();
    }

    // Return the results as a tuple
    return {time_new, time_old, time_baseline};
}

// Function to handle switch-case for MPFR precision without known dataSize
std::tuple<double, double, double> handle_mpfr_precision(int precision, const Eigen::MatrixXd &ptsA, 
                                                         const Eigen::MatrixXd &ptsB, 
                                                         const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    switch (precision) {
        case 16: return run_mpfr_benchmark<16>(ptsA, ptsB, lattitudes, X, Y);
        case 17: return run_mpfr_benchmark<17>(ptsA, ptsB, lattitudes, X, Y);
        case 18: return run_mpfr_benchmark<18>(ptsA, ptsB, lattitudes, X, Y);
        case 19: return run_mpfr_benchmark<19>(ptsA, ptsB, lattitudes, X, Y);
        case 20: return run_mpfr_benchmark<20>(ptsA, ptsB, lattitudes, X, Y);
        case 21: return run_mpfr_benchmark<21>(ptsA, ptsB, lattitudes, X, Y);
        case 22: return run_mpfr_benchmark<22>(ptsA, ptsB, lattitudes, X, Y);
        case 23: return run_mpfr_benchmark<23>(ptsA, ptsB, lattitudes, X, Y);
        case 24: return run_mpfr_benchmark<24>(ptsA, ptsB, lattitudes, X, Y);
        case 25: return run_mpfr_benchmark<25>(ptsA, ptsB, lattitudes, X, Y);
        case 26: return run_mpfr_benchmark<26>(ptsA, ptsB, lattitudes, X, Y);
        case 27: return run_mpfr_benchmark<27>(ptsA, ptsB, lattitudes, X, Y);
        case 28: return run_mpfr_benchmark<28>(ptsA, ptsB, lattitudes, X, Y);
        case 29: return run_mpfr_benchmark<29>(ptsA, ptsB, lattitudes, X, Y);
        case 30: return run_mpfr_benchmark<30>(ptsA, ptsB, lattitudes, X, Y);
        case 31: return run_mpfr_benchmark<31>(ptsA, ptsB, lattitudes, X, Y);
        case 32: return run_mpfr_benchmark<32>(ptsA, ptsB, lattitudes, X, Y);
        default:
            std::cerr << "Unsupported MPFR precision: " << precision << std::endl;
            return {0.0, 0.0, 0.0}; // Return a default tuple in case of an error
    }
}


// Helper function to run Quadruple benchmark based on precision
std::tuple<double, double, double> run_quadruple_benchmark(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    // Get the intersection functions for the three different equations (new, old, baseline)

    // To store the timings for each method
    double time_new = 0.0, time_old = 0.0, time_baseline = 0.0;

    // Measure time for the "new" method
    {
        auto start_new = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < numTests; ++i) {
            // Do not use std::move on X and Y here since we need to reuse them in the loop
            compute_values_new_cast<mppp::real128>(ptsA, ptsB, lattitudes, X, Y);
        }
        std::chrono::duration<double> elapsed_new = std::chrono::high_resolution_clock::now() - start_new;
        time_new = elapsed_new.count();
    }

    // Measure time for the "old" method
    {
        auto start_old = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < numTests; ++i) {
            compute_values_old_cast<mppp::real128>(ptsA, ptsB, lattitudes, X, Y);
        }
        std::chrono::duration<double> elapsed_old = std::chrono::high_resolution_clock::now() - start_old;
        time_old = elapsed_old.count();
    }

    // Measure time for the "baseline" method
    {
        auto start_baseline = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < numTests; ++i) {
            compute_values_baseline_cast<mppp::real128>(ptsA, ptsB, lattitudes, X, Y);
        }
        std::chrono::duration<double> elapsed_baseline = std::chrono::high_resolution_clock::now() - start_baseline;
        time_baseline = elapsed_baseline.count();
    }

    // Return the results as a tuple
    return {time_new, time_old, time_baseline};
}


#endif // SCALAR_TYPE_BENCHMARK_HH