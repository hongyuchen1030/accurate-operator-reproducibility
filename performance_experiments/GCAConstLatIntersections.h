//
// Created by hongyu chen on 4/8/24.
//

#ifndef CODES_GCACONSTLATINTERSECTIONS_VECTORIZED_H
#define CODES_GCACONSTLATINTERSECTIONS_VECTORIZED_H

#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <Eigen/Dense>
// #include <boost/multiprecision/mpfr.hpp>
// #include <boost/multiprecision/gmp.hpp>
#include <boost/math/tools/roots.hpp>
#include <type_traits>

template<typename T, size_t N>
using VN_T = Eigen::Matrix<T, N, 1>;
template<typename T>
using V3_T = VN_T<T, 3>;


template<typename T>
static V3_T<T> simd_cross(const V3_T<T>& v1, const V3_T<T>& v2){    
    V3_T<T> n;
    n(0) = v1.coeff(1) * v2.coeff(2) - v1.coeff(2) * v2.coeff(1);  // x = Ay * Bz - Az * By
    n(1) = v1.coeff(2) * v2.coeff(0) - v1.coeff(0) * v2.coeff(2);  // y = Az * Bx - Ax * Bz
    n(2) = v1.coeff(0) * v2.coeff(1) - v1.coeff(1) * v2.coeff(0);  // z = Ax * By - Ay * Bx
    return n;
}


template<typename T>
static std::tuple<T, T> gca_constLat_intersection_coordinates_baselineEqn_performance(const V3_T<T>& pointA,const  V3_T<T>& pointB, const T &constZ) {

    // Calculate the cross product, which gives us the vector n and normalize it
    V3_T<T> n = simd_cross(pointA, pointB);

    //Normalize n
    n.normalize();

    // Extract nx, ny, nz components from vector n
    T nx = n[0];
    T ny = n[1];
    T nz = n[2];


    T nx_squared = nx * nx;
    T ny_squared = ny * ny;
    T nz_squared = nz * nz;
    T nx_squared_plus_ny_squared = nx_squared + ny_squared;

    // Calculate s (as 's' in the formula)
    T s = sqrt(nx_squared_plus_ny_squared - constZ * constZ);

    // Calculate p_x and p_y, which are the x and y coordinates of the intersection
    T p_x = - (constZ * nx * nz + s * ny) / 1.0 - nz_squared;
    T p_y = - (constZ * ny * nz - s * nx) / 1.0 - nz_squared; // Only consider the minus


    return {p_x, p_y};
}


template<typename T>
static std::tuple<T, T> gca_constLat_intersection_coordinates_oldEqn_performance(const V3_T<T>& pointA, const V3_T<T>& pointB, const T &constZ) {

    // Calculate the cross product, which gives us the vector n and normalize it
    V3_T<T> n = simd_cross(pointA, pointB);
    n.normalize();

    // Extract nx, ny, nz components from vector n
    T nx = n[0];
    T ny = n[1];
    T nz = n[2];


    T nx_squared = nx * nx;
    T ny_squared = ny * ny;
    T nz_squared = nz * nz;
    T nx_squared_plus_ny_squared = nx_squared + ny_squared;




    // Calculate s (as 's' in the formula)
    T s = sqrt(nx_squared_plus_ny_squared - constZ * constZ);

    // Calculate p_x and p_y, which are the x and y coordinates of the intersection
    T p_x = - (constZ * nx * nz + s * ny) / nx_squared_plus_ny_squared;
    T p_y = - (constZ * ny * nz - s * nx) / nx_squared_plus_ny_squared; // Only consider the minus


    return {p_x, p_y};
}

template<typename T>
static std::tuple<T, T> gca_constLat_intersection_coordinates_newEqn_performance(const V3_T<T>& pointA, const V3_T<T>& pointB, const T &constZ) {
    // Calculate the cross product, which gives us the vector n
    V3_T<T> n = simd_cross(pointA, pointB);

    // Extract nx, ny, nz components from vector n
    T nx = n[0];
    T ny = n[1];
    T nz = n[2];



    // Calculate s (as 's' in the formula)
    T nx_squared = nx * nx;
    T ny_squared = ny * ny;
    T nz_squared = nz * nz;
    T nx_squared_plus_ny_squared = nx_squared + ny_squared;
    T norm_n_squared = nx_squared_plus_ny_squared + nz_squared; // ||n||^2


    T s_tilde = sqrt(nx_squared_plus_ny_squared - norm_n_squared * constZ * constZ);

    // T s_tilde = sqrt(nx_squared_plus_ny_squared - constZ * constZ);

    // Calculate p_x and p_y, which are the x and y coordinates of the intersection
    T p_x = - (constZ * nx * nz + s_tilde * ny)/nx_squared_plus_ny_squared;

    T p_y = - (constZ * ny * nz - s_tilde * nx)/nx_squared_plus_ny_squared; // Only consider the minus



    return {p_x, p_y};
}



class AccurateIntersection {
public:
    template <typename T>
    static std::tuple<T, T> gca_constLat_intersection_accurate_performance(const std::array<T, 3>& x1, const std::array<T, 3>& x2, T z0) {
        // Calculate nx, ny, nz and their errors
        auto [nx, enx] = fmms_re(x1[1], x2[2], x1[2], x2[1]); // Component x
        auto [ny, eny] = fmms_re(x1[2], x2[0], x1[0], x2[2]); // Component y
        auto [nz, enz] = fmms_re(x1[0], x2[1], x1[1], x2[0]); // Component z

        // Calculate sum of squares for nx and ny
        std::array<T, 2> nx_ny = {nx, ny};
        std::array<T, 2> enx_eny = {enx, eny};
        auto [nx_squared_plus_ny_squared_pre, nx_squared_plus_ny_squared_err_pre] = sum_of_squares_re(nx_ny, enx_eny);

        auto [nx_squared_plus_ny_squared, nx_squared_plus_ny_squared_err] = two_sum(nx_squared_plus_ny_squared_pre,  nx_squared_plus_ny_squared_err_pre);

        // Sum of squares for nx, ny, and nz
        std::array<T, 3> nxyz = {nx, ny, nz};
        std::array<T, 3> enxyz = {enx, eny, enz};
        auto [norm_n_squared_pre, norm_n_squared_err_pre] = sum_of_squares_re(nxyz, enxyz);

        auto [norm_n_squared, norm_n_squared_err] = two_sum(norm_n_squared_pre, norm_n_squared_err_pre);
        // TwoProd for z0
        auto [C, c] = two_prod_fma(z0, z0);

        // Compute dot product for J with itself and C
        std::array<T, 4> norm_n_squared_array = {norm_n_squared, norm_n_squared, norm_n_squared_err, norm_n_squared_err};
        std::array<T, 4> C_array = {C, c, C, c};
        auto [J, j] = dot_fma_re(norm_n_squared_array, C_array);

        // TwoSum for -J and nx_squared_plus_ny_squared
        auto [s2, f] = two_sum(T(-J), nx_squared_plus_ny_squared);

        // VecSum for error correction
        T es2 = -j + nx_squared_plus_ny_squared_err + f;

        // Apply two_sum to combine s2 and es2
        auto [s2_, es2_] = two_sum(s2, es2);

        // Accurate sqrt for s2 and es2
        auto [s, es] = acc_sqrt_re(s2_, es2_);

        // Calculate terms for p_x
        auto [A_px, ea1_px] = two_prod_fma(nx, nz);
        T ea2_px = nx * enz;
        T ea3_px = nz * enx;
        T ea_px = ea1_px + ea2_px + ea3_px;

        std::array<T, 6> vec1x = {A_px, ea_px, ny, eny, ny, eny};
        std::array<T, 6> vec2x = {z0, z0, s, s, es, es};
        T D = dot_fma(vec1x, vec2x);

        // Calculate terms for p_y
        auto [H, eh1] = two_prod_fma(ny, nz);
        T eh2 = ny * enz;
        T eh3 = nz * eny;
        T eh = eh1 + eh2 + eh3;

        std::array<T, 6> vec1y = {H, eh, -nx, -enx, -nx, -enx};
        std::array<T, 6> vec2y = {z0, z0, s, s, es, es};

        T E = dot_fma(vec1y, vec2y);

        T res_x = -D/nx_squared_plus_ny_squared;
        T res_y = -E/nx_squared_plus_ny_squared;
        return {res_x, res_y};
    }

private:

    // Helper function for two_prod_fma using template T
    template <typename T>
    static std::tuple<T, T> two_prod_fma(T a, T b) {
        T x = a * b;
        T neg_x = -x;
        T y = simd_fma(a, b, neg_x);
        return std::make_tuple(x, y);
    }

    template<typename T>
    static std::tuple<T, T> two_sum(const T& a, const T& b) {
        T x = a + b;
        T z = x - a;
        return std::make_tuple(x, (a - (x - z)) + (b - z));
    }

    // FMA-based product difference using template T
    template <typename T>
    static std::tuple<T, T> fmms_re(T a, T b, T c, T d) {
        auto [p1, s1] = two_prod_fma(a, b);  
        T neg_d = -d;
        auto [h2, r2] = two_prod_fma(c, neg_d);
        auto [p2, q2] = two_sum(p1, h2);    
        T s2 = s1 + q2 + r2;
       
        return  std::make_tuple(p2, s2);
    }

    // Template function for two_square_fma with generic type T
    template <typename T>
    static std::tuple<T, T> two_square_fma(T a) {
        T x = a * a;
        T neg_x = -x;
        T y = simd_fma(a, a, neg_x);
        return std::make_tuple(x, y);
    }

    template <typename T>
    static std::tuple<T, T> fast_two_sum(T a, T b) {
        T x = a + b;
        T y = a - x + b;
        return  std::make_tuple(x, y);
    }

    // Accurate sqrt calculation with error compensation using template T
    template <typename T>
    static std::tuple<T, T> acc_sqrt_re(T T_val, T t) {
        T P = sqrt(T_val);
        auto [H, h] = two_square_fma(P);
        T r =  (T_val - H - h) + t;
        T p = 0.5 * r * P;

        return std::make_tuple(P, p);
    }

    template <typename T, std::size_t N>
    static std::tuple<T, T> dot_fma_re(const std::array<T, N>& v1, const std::array<T, N>& v2) {
        auto [s, c] = two_prod_fma(v1[0], v2[0]);
        for (std::size_t i = 1; i < N; ++i) {
            auto [p, pi] = two_prod_fma(v1[i], v2[i]);
            T sigma;
            std::tie(s, sigma) = two_sum(s, p); // Potential opportunity for speedup: do we really need this accuracy in the error term?
            c += pi + sigma;
        }
        return std::make_tuple(s, c);  // Return the dot product and compensation term
    }

    template <typename T, std::size_t N>
    static T dot_fma(const std::array<T, N>& v1, const std::array<T, N>& v2) {
        auto [s, c] = dot_fma_re(v1, v2);
        return s + c;
    }

    template <typename T>
    static std::tuple<T, T> sum_non_neg(T A_high, T A_low, T B_high, T B_low) {
        // High part of the sum using two_sum
        auto [H, h] = two_sum(A_high, B_high);

        // Low part of the sum
        T c = A_low + B_low;
        T d = h + c;

        // Apply fast_two_sum for the final summation
        return fast_two_sum(H, d);
    }

    template <typename T, std::size_t N>
    static std::tuple<T, T> sum_of_squares_re(const std::array<T, N>& vals, const std::array<T, N>& errs) {
        // Initialize S and s based on type T
        T S;
        T s;
        S = 0.0;
        s = 0.0;

        for (std::size_t i = 0; i < N; ++i) {
            auto [P, p] = two_prod_fma(vals[i], vals[i]);  // Compute product and error for vals[i]^2
            std::tie(S, s) = sum_non_neg(S, s, P, p);  // Sum with compensation
        }

        // Compute the dot product between the values and errors
        T R = dot_fma(vals, errs);

        return std::make_tuple(S, T(2 * R + s)); // Note: 2 * R is exact assuming no overflow
    }
};
#endif //CODES_GCACONSTLATINTERSECTIONS_H
