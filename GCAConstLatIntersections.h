//
// Created by hongyu chen on 4/8/24.
//

#ifndef CODES_GCACONSTLATINTERSECTIONS_H
#define CODES_GCACONSTLATINTERSECTIONS_H

#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <Eigen/Dense>
#include <boost/math/tools/roots.hpp>
#include <type_traits>
#include "simd_fma.hh"
template<typename T, size_t N>
using VN_T = Eigen::Matrix<T, N, 1>;
template<typename T>
using V3_T = VN_T<T, 3>;



// Helper function to initialize T to zero
template <typename T>
T initialize_zero(const T& example) {
    if constexpr (std::is_scalar_v<T>) {
        return 0.0;  // Return scalar zero
    } else {
        return T::Zero(example.size());  // Return zero vector (for Eigen::VectorXd)
    }
}


template<typename T>
static V3_T<T> simd_cross(const V3_T<T>& v1, const V3_T<T>& v2){    
    V3_T<T> n;
    n(0) = v1.coeff(1) * v2.coeff(2) - v1.coeff(2) * v2.coeff(1);  // x = Ay * Bz - Az * By
    n(1) = v1.coeff(2) * v2.coeff(0) - v1.coeff(0) * v2.coeff(2);  // y = Az * Bx - Ax * Bz
    n(2) = v1.coeff(0) * v2.coeff(1) - v1.coeff(1) * v2.coeff(0);  // z = Ax * By - Ay * Bx
    return n;
}

template<typename T>
struct intermediate_values_compare {
    T nx;
    T nNorm;
    T sSquare;
    T s;
    T znxnz;
    T sny;
    T znxnz_sny;
    T nxSquarenySquare;
    T px;

    // Default constructor
    intermediate_values_compare()
            : nx(T()), nNorm(T()), sSquare(T()), s(T()), znxnz(T()), sny(T()), znxnz_sny(T()), nxSquarenySquare(T()), px(T()) {}

    // Parameterized constructor
    intermediate_values_compare(T nx, T nNorm, T sSquare, T s, T znxnz, T sny, T znxnz_sny, T nxSquarenySquare, T px)
            : nx(nx), nNorm(nNorm), sSquare(sSquare), s(s), znxnz(znxnz), sny(sny), znxnz_sny(znxnz_sny), nxSquarenySquare(nxSquarenySquare), px(px) {}
};




struct intermediate_values_our {
    std::tuple<double, double> nx_pair;
    std::tuple<double, double> nNorm_pair; // Changed from nNorm_quad to nNorm_pair
    std::tuple<double, double> sSquare_pair;
    std::tuple<double, double> s_pair;
    std::tuple<double, double> znxnz_pair;
    std::tuple<double, double> sny_pair;
    double znxnz_sny;
    double nxSquarenySquare;
    double px;

    // Default constructor
    intermediate_values_our()
        : nx_pair(0.0, 0.0), nNorm_pair(0.0, 0.0), sSquare_pair(0.0, 0.0),
          s_pair(0.0, 0.0), znxnz_pair(0.0, 0.0), sny_pair(0.0, 0.0),
          znxnz_sny(0.0), nxSquarenySquare(0.0), px(0.0) {}

    // Parameterized constructor
    intermediate_values_our(std::tuple<double, double> nx_pair,
                            std::tuple<double, double> nNorm_pair, // Changed from nNorm_quad to nNorm_pair
                            std::tuple<double, double> sSquare_pair,
                            std::tuple<double, double> s_pair,
                            std::tuple<double, double> znxnz_pair,
                            std::tuple<double, double> sny_pair,
                            double znxnz_sny,
                            double nxSquarenySquare,
                            double px)
        : nx_pair(nx_pair), nNorm_pair(nNorm_pair), sSquare_pair(sSquare_pair), s_pair(s_pair),
          znxnz_pair(znxnz_pair), sny_pair(sny_pair), znxnz_sny(znxnz_sny), nxSquarenySquare(nxSquarenySquare), px(px) {}
};



template<typename T>
static std::tuple<T, T> gca_constLat_intersection_coordinates_baselineEqn(V3_T<T>& pointA, V3_T<T>& pointB, T constZ) {
    // Ensure |z1| >= |z2| by swapping if necessary
//    if (abs(pointA[2]) < std::abs(pointB[2])) {
//        pointA.swap(pointB);
//    }
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

    // Calculate s (as 's' in the formula)
    T s = sqrt(nx_squared + ny_squared - constZ * constZ);
    T denominator = 1.0 - nz_squared;//1.0 - nz_squared;

    // Calculate p_x and p_y, which are the x and y coordinates of the intersection
    T p_x = - (constZ * nx * nz + s * ny) / denominator;
    T p_y = - (constZ * ny * nz - s * nx) / denominator; // Only consider the minus


    return {p_x, p_y};
}

template<typename T>
static std::tuple<T, T> gca_constLat_intersection_coordinates_baselineEqn_performance(V3_T<T>& pointA, V3_T<T>& pointB, T constZ) {
    // Ensure |z1| >= |z2| by swapping if necessary
//    if (abs(pointA[2]) < std::abs(pointB[2])) {
//        pointA.swap(pointB);
//    }
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
static std::tuple<T, T> gca_constLat_intersection_coordinates_oldEqn(V3_T<T>& pointA, V3_T<T>& pointB, T constZ) {
    // Ensure |z1| >= |z2| by swapping if necessary
//    if (abs(pointA[2]) < std::abs(pointB[2])) {
//        pointA.swap(pointB);
//    }
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
    // T denominator = nx_squared+ny_squared;//1.0 - nz_squared;

    // Calculate p_x and p_y, which are the x and y coordinates of the intersection
    T p_x = - (constZ * nx * nz + s * ny) / nx_squared_plus_ny_squared;
    T p_y = - (constZ * ny * nz - s * nx) / nx_squared_plus_ny_squared; // Only consider the minus


    return {p_x, p_y};
}


template<typename T>
static std::tuple<T, T> gca_constLat_intersection_coordinates_newEqn(
        const V3_T<T>& pointA, const V3_T<T>& pointB, T constZ) {
    intermediate_values_compare<T> dummy;
    return gca_constLat_intersection_coordinates_newEqn(pointA, pointB, constZ, dummy);
}


template<typename T>
static std::tuple<T, T> gca_constLat_intersection_coordinates_newEqn(const V3_T<T>& pointA, const V3_T<T>& pointB, T constZ,
                                                                     intermediate_values_compare<T>& intermediateValuesCompare) {
    // Calculate the cross product, which gives us the vector n
    V3_T<T> n = simd_cross(pointA, pointB);

    //Normalize n
    n.normalize();


    // Extract nx, ny, nz components from vector n
    T nx = n[0];
    T ny = n[1];
    T nz = n[2];

    //print out the calculated nx, ny, nz
    intermediateValuesCompare.nx = nx;


    // Calculate s (as 's' in the formula)
    T nx_squared = nx * nx;
    T ny_squared = ny * ny;
    T nz_squared = nz * nz;
    T nx_squared_plus_ny_squared = nx_squared + ny_squared;

    T norm_n_squared = nx_squared_plus_ny_squared + nz_squared; // ||n||^2

    T s_tilde = sqrt(nx_squared_plus_ny_squared - norm_n_squared * constZ * constZ);

    // print out the calculated norm_n_squared
    intermediateValuesCompare.nxSquarenySquare = nx_squared+ ny_squared;
    intermediateValuesCompare.nNorm = norm_n_squared;
    intermediateValuesCompare.sSquare = (nx_squared_plus_ny_squared - norm_n_squared * constZ * constZ);
    intermediateValuesCompare.s = s_tilde;
    intermediateValuesCompare.znxnz_sny = (constZ * nx * nz + s_tilde * ny);

    // Calculate p_x and p_y, which are the x and y coordinates of the intersection
    T p_x = - (constZ * nx * nz + s_tilde * ny)/nx_squared_plus_ny_squared;

    intermediateValuesCompare.znxnz = (constZ * nx * nz);

    intermediateValuesCompare.sny = (s_tilde * ny);
    T p_y = - (constZ * ny * nz - s_tilde * nx)/nx_squared_plus_ny_squared; // Only consider the minus
    intermediateValuesCompare.px = p_x;

    return {p_x, p_y};
}



class AccurateIntersection {
public:
    // Constructor (if needed)
    AccurateIntersection() = default;

    std::tuple<double, double> gca_constLat_intersection_accurate(const std::array<double, 3>& x1, const std::array<double, 3>& x2, double z0){
        intermediate_values_our intermediateValuesOur;
        return gca_constLat_intersection_accurate(x1, x2, z0, intermediateValuesOur);
    }

    std::tuple<double, double> gca_constLat_intersection_accurate(const std::array<double, 3>& x1, const std::array<double, 3>& x2, double z0, intermediate_values_our& intermediateValuesOur) {
        // Step 1: Calculate nx, ny, nz and their errors using the fmms_re helper
        auto [nx, enx] = fmms_re(x1[1], x2[2], x1[2], x2[1]);  // Component x
        auto [ny, eny] = fmms_re(x1[2], x2[0], x1[0], x2[2]);  // Component y
        auto [nz, enz] = fmms_re(x1[0], x2[1], x1[1], x2[0]);  // Component z

        intermediateValuesOur.nx_pair = {nx, enx};

        // Step 2: Calculate sum of squares for nx and ny
        auto [nx_squared_plus_ny_squared_pre, nx_squared_plus_ny_squared_err_pre] = sum_of_squares_re(std::array<double, 2>{nx, ny}, std::array<double, 2>{enx, eny});
        auto [nx_squared_plus_ny_squared, nx_squared_plus_ny_squared_err] = two_sum(nx_squared_plus_ny_squared_pre, nx_squared_plus_ny_squared_err_pre);
        intermediateValuesOur.nxSquarenySquare = nx_squared_plus_ny_squared;

        // Step 3: Calculate sum of squares for the full vector including nz
        auto [norm_n_squared_pre, norm_n_squared_err_pre] = sum_of_squares_re(std::array<double, 3>{nx, ny, nz}, std::array<double, 3>{enx, eny, enz});
        auto [norm_n_squared, norm_n_squared_err] = two_sum(norm_n_squared_pre, norm_n_squared_err_pre);
        intermediateValuesOur.nNorm_pair = {norm_n_squared, norm_n_squared_err};

        // Step 4: TwoProd for z0
        auto [C, c] = two_prod_fma(z0, z0);

        // Step 5: Calculate J and j using dot_fma_re with norm_n_squared and C
        auto [J, j] = dot_fma_re(std::array<double, 4>{norm_n_squared, norm_n_squared, norm_n_squared_err, norm_n_squared_err}, std::array<double, 4>{C, c, C, c});

        // Step 6: TwoSum for -J and nx_squared_plus_ny_squared
        auto [s2, f] = two_sum(-J, nx_squared_plus_ny_squared);

        // Step 7: Vector sum for -j, nx_squared_plus_ny_squared_err, and f
        double es2 = vec_sum(std::array<double, 3>{-j, nx_squared_plus_ny_squared_err, f});

        // Step 8: Apply two_sum to combine s2 and es2
        auto [s2_, es2_] = two_sum(s2, es2);
        intermediateValuesOur.sSquare_pair = {s2_, es2_};

        // Step 9: Calculate sqrt for s2 and es2 using acc_sqrt_re
        auto [s, es] = acc_sqrt_re(s2_, es2_);
        intermediateValuesOur.s_pair = {s, es};

        // Step 10: Calculate the terms for s * nx * nz using two_prod_fma
        auto [A_px, ea1_px] = two_prod_fma(nx, nz);
        double ea2_px = nx * enz;
        double ea3_px = nz * enx;
        double ea_px = ea1_px + ea2_px + ea3_px;

        std::array<double, 6> vec1x = {A_px, ea_px, ny, eny, ny, eny};
        std::array<double, 6> vec2x = {z0, z0, s, s, es, es};

        // Step 11: Use dot_fma to calculate D
        double D = dot_fma(vec1x, vec2x);

        // Step 12: Store intermediate values
        intermediateValuesOur.znxnz_pair = dot_fma_re(std::array<double, 2>{A_px, ea_px}, std::array<double, 2>{z0, z0});
        intermediateValuesOur.sny_pair = dot_fma_re(std::array<double, 4>{ny, eny, ny, eny}, std::array<double, 4>{s, s, es, es});
        intermediateValuesOur.znxnz_sny = D;

        // Step 13: Calculate the terms for p_y
        auto [H, eh1] = two_prod_fma(ny, nz);
        double eh2 = ny * enz;
        double eh3 = nz * eny;
        double eh = eh1 + eh2 + eh3;

        std::array<double, 6> vec1y = {H, eh, -nx, -enx, -nx, -enx};
        std::array<double, 6> vec2y = {z0, z0, s, s, es, es};

        // Step 14: Use dot_fma to calculate E
        double E = dot_fma(vec1y, vec2y);

        // Step 15: Calculate final results res_x and res_y
        double res_x = -D / nx_squared_plus_ny_squared;
        intermediateValuesOur.px = res_x;
        double res_y = -E / nx_squared_plus_ny_squared;

        return {res_x, res_y};
    }


private:
    // Helper function for two_prod_fma using template T
    template <typename T>
    std::tuple<T, T> two_prod_fma(T a, T b) {
        T x = a * b;
        T neg_x = -x;
        T y = simd_fma(a, b, neg_x);
        return std::make_tuple(x, y);
    }

    // Helper function for two_sum using template T
    template <typename T>
    std::tuple<T, T> two_sum(const T& a, const T& b) {
        T x = a + b;
        T z = x - a;
        T y = (a - (x - z)) + (b - z);
        return std::make_tuple(x, y);
    }

    // FMA-based product difference using template T
    template <typename T>
    std::tuple<T, T> fmms_re(T a, T b, T c, T d) {
        auto [p1, s1] = two_prod_fma(a, b);  
        T neg_d = -d;
        auto [h2, r2] = two_prod_fma(c, neg_d);
        auto [p2, q2] = two_sum(p1, h2);    
        T s2 = s1 + (q2 + r2);
       
        return  std::make_tuple(p2, s2);
    }

    // Helper function to calculate (a * b) - (c * d) using FMA
    template <typename T>
    T fmms(T a, T b, T c, T d) {
        T cd = c * d;
        T neg_c = -c;
        T err = simd_fma(neg_c, d, cd);
        T neg_cd = -cd;
        T dop = simd_fma(a, b, neg_cd);
        return dop + err;
    }

    // Calculate the cross product of two 3D vectors utilizing the FMA operation
    template <typename T>
    std::array<T, 3> cross_fma(const std::array<T, 3>& v1, const std::array<T, 3>& v2) {
        std::array<T, 3> result;
        result[0] = fmms(v1[1], v2[2], v1[2], v2[1]);
        result[1] = fmms(v1[2], v2[0], v1[0], v2[2]);
        result[2] = fmms(v1[0], v2[1], v1[1], v2[0]);
        return result;
    }

    // Template function for two_square_fma with generic type T
    template <typename T>
    std::tuple<T, T> two_square_fma(T a) {
        T x = a * a;
        T neg_x = -x;
        T y = simd_fma(a, a, neg_x);
        return std::make_tuple(x, y);
    }

    // Fast two-sum helper function using template T
    template <typename T>
    std::tuple<T, T> fast_two_sum(T a, T b) {
        T x = a + b;
        T y = (a - x) + b;
       
        return  std::make_tuple(x, y);
    }

    // Accurate sqrt calculation with error compensation using template T
    template <typename T>
    std::tuple<T, T> acc_sqrt_re(T T_val, T t) {
        T P = sqrt(T_val);
        T H, h;
        std::tie(H, h) = two_square_fma(P);
        T r = (T_val - H) - h;
        r += t;
        T p = r / (2.0 * P);
        
        return std::make_tuple(P, p);
    }

    template <typename T, std::size_t N>
    T dot_fma(const std::array<T, N>& v1, const std::array<T, N>& v2) {
        T s, c;
        std::tie(s, c) = two_prod_fma(v1[0], v2[0]);
        for (std::size_t i = 1; i < N; ++i) {
            T p, pi;
            std::tie(p, pi) = two_prod_fma(v1[i], v2[i]);
            T sigma;
            std::tie(s, sigma) = two_sum(s, p);
            c += pi + sigma;
        }
        return s + c;
}



    // General template for dot_fma_re when T is a scalar (double)
    template <typename T, std::size_t N>
    std::tuple<T, T> dot_fma_re(const std::array<T, N>& v1, const std::array<T, N>& v2) {
        T s, c;
        std::tie(s, c) = two_prod_fma(v1[0], v2[0]);
        for (std::size_t i = 1; i < N; ++i) {
            T p, pi;
            std::tie(p, pi) = two_prod_fma(v1[i], v2[i]);
            T sigma;
            std::tie(s, sigma) = two_sum(s, p);
            c += pi + sigma;
        }
        return std::make_tuple(s, c);  // Return the dot product and compensation term
    }


    // General template for sum_non_neg
    template <typename T>
    std::tuple<T, T> sum_non_neg(const std::tuple<T, T>& A, const std::tuple<T, T>& B) {
        T A_high, A_low, B_high, B_low;
        std::tie(A_high, A_low) = A;
        std::tie(B_high, B_low) = B;

        // High part of the sum using two_sum
        auto [H, h] = two_sum(A_high, B_high);

        // Low part of the sum
        T c = A_low + B_low;
        T d = h + c;

        // Apply fast_two_sum for the final summation
        auto [S, s] = fast_two_sum(H, d);

        return std::make_tuple(S, s);
    }




    // General template for vec_sum when T is double
    template <typename T, std::size_t N>
    T vec_sum(const std::array<T, N>& arr) {
        T sum;  // Initialize sum with the first element
        sum = arr[0];
        T compensator; // Initialize compensator to zero
        compensator = 0.0;

        for (std::size_t i = 1; i < N; ++i) {
            auto [partial_sum, err] = two_sum(sum, arr[i]);  // Perform two_sum to maintain accuracy
            sum = partial_sum;
            compensator += err;  // Accumulate the error term
        }

        return sum + compensator;  // Return the compensated sum
    }


    // General template for sum_of_squares_re when T is double
    template <typename T, std::size_t N>
    std::tuple<T, T> sum_of_squares_re(const std::array<T, N>& vals, const std::array<T, N>& errs) {
        // Initialize S and s based on type T
        T S;
        T s;
        S = 0.0;
        s = 0.0;

        for (std::size_t i = 0; i < N; ++i) {
            auto [P, p] = two_prod_fma(vals[i], vals[i]);  // Compute product and error for vals[i]^2
            auto [new_S, new_s] = sum_non_neg(std::make_tuple(S, s), std::make_tuple(P, p));  // Sum with compensation
            S = new_S;
            s = new_s;
        }

        // Compute the dot product between the values and errors
        T R = dot_fma(vals, errs);

        // Error term
        T err = 2 * R + s;

        return std::make_tuple(S, err);  // Return the sum of squares and the error term
    }



};


#endif //CODES_GCACONSTLATINTERSECTIONS_H
