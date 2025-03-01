#ifndef VECTORIZED_BENCHMARK_TEMPLATE_HH
#define VECTORIZED_BENCHMARK_TEMPLATE_HH

#include <Eigen/Dense>

// Include other headers required for your code
#include "simd_fma.hh"
#include "GCAConstLatIntersections.h"
#include <tbb/parallel_for.h>
template<typename T, size_t N>
using VN_T = Eigen::Matrix<T, N, 1>;
template<typename T>
using V3_T = VN_T<T, 3>;

template <typename T>
constexpr size_t vecWidth() {
    if constexpr (std::is_arithmetic_v<T>) return 1;
    else return T::RowsAtCompileTime;
}

template<size_t VEC_WIDTH>
using Vec = Eigen::Array<double, VEC_WIDTH, 1>;

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
void compute_values(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    const size_t size = ptsA.rows();
    constexpr size_t stride = vecWidth<T>();
    assert((size % stride == 0) && "Size must be a multiple of the vector width");

    if ((size != ptsB.rows()) || (size != lattitudes.rows()) || (size != X.rows()) || (size != Y.rows())) {
        throw std::runtime_error("All passed arrays must have the same size");
    }

    auto processChunk = [&](size_t chunk_i) {
        const size_t i = chunk_i * stride;
        V3_T<T> x1, x2;
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

        auto [px, py] = gca_constLat_intersection_coordinates_newEqn_performance(x1, x2, z0);

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
                        processChunk(i);
                });
    } else {
        for (size_t i = 0; i < size / stride; ++i) { processChunk(i); }
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



template <typename T, bool Parallel = false>
double run_benchmark(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    auto start = std::chrono::high_resolution_clock::now();

    constexpr size_t numTests = 100;
    for (size_t i = 0; i < numTests; ++i) {
         compute_values<T, Parallel>(ptsA, ptsB, lattitudes, X, Y);
    }

    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;

    return elapsed.count();
}

#endif  // VECTORIZED_BENCHMARK_TEMPLATE_HH
