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


    n.normalize();


    T nx = n[0];
    T ny = n[1];
    T nz = n[2];


    T nx_squared = nx * nx;
    T ny_squared = ny * ny;
    T nz_squared = nz * nz;
    T denominator = 1.0 - nz_squared;

    // Calculate s (as 's' in the formula)
    T s = sqrt(denominator - constZ * constZ);

    // Calculate p_x and p_y, which are the x and y coordinates of the intersection
    T p_x = - (constZ * nx * nz + s * ny) / denominator;
    T p_y = - (constZ * ny * nz - s * nx) / denominator; // Only consider the minus


    return {p_x, p_y};
}


template<typename T>
static std::tuple<T, T> gca_constLat_intersection_coordinates_oldEqn_performance(const V3_T<T>& pointA, const V3_T<T>& pointB, const T &constZ) {

    // Calculate the cross product, which gives us the vector n and normalize it
    V3_T<T> n = simd_cross(pointA, pointB);
    n.normalize();


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

    V3_T<T> n = simd_cross(pointA, pointB);


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


    // Calculate p_x and p_y, which are the x and y coordinates of the intersection
    T p_x = - (constZ * nx * nz + s_tilde * ny)/nx_squared_plus_ny_squared;

    T p_y = - (constZ * ny * nz - s_tilde * nx)/nx_squared_plus_ny_squared; // Only consider the minus



    return {p_x, p_y};
}



class AccurateIntersection {
public:
    template <typename T>
    static std::tuple<T, T> gca_constLat_intersection_accurate_performance(const std::array<T, 3>& x1, const std::array<T, 3>& x2, T z0) {

        auto [nx, enx] = AccuDOP(x1[1], x2[2], x1[2], x2[1]); 
        auto [ny, eny] = AccuDOP(x1[2], x2[0], x1[0], x2[2]); 
        auto [nz, enz] = AccuDOP(x1[0], x2[1], x1[1], x2[0]); 


        std::array<T, 2> nx_ny = {nx, ny};
        std::array<T, 2> enx_eny = {enx, eny};
        auto [S_2, s_2_err] = SumOfSquaresC(nx_ny, enx_eny);


        std::array<T, 3> nxyz = {nx, ny, nz};
        std::array<T, 3> enxyz = {enx, eny, enz};
        auto [S_3, s_3_err] = SumOfSquaresC(nxyz, enxyz);
        

        auto [C, c] = two_prod_fma(z0, z0);


        std::array<T, 4> norm_n_squared_array = {S_3, S_3, s_3_err, s_3_err};
        std::array<T, 4> C_array = {C, c, C, c};
        auto [D, e_D] = CompDotC(norm_n_squared_array, C_array);


        auto [E, e_err] = two_sum(T(-D), S_2);


        T es2 = -e_D + s_2_err + e_err;


        auto [s2_, es2_] = two_sum(E, es2);


        auto [s, es] = acc_sqrt_re(s2_, es2_);


        auto [A_px, ea1_px] = two_prod_fma(nx, nz);
        T ea2_px = nx * enz;
        T ea3_px = nz * enx;
        T ea_px = ea1_px + ea2_px + ea3_px;

        std::array<T, 6> vec1x = {A_px, ea_px, ny, eny, ny, eny};
        std::array<T, 6> vec2x = {z0, z0, s, s, es, es};
        T x = dot_fma(vec1x, vec2x);


        auto [H, eh1] = two_prod_fma(ny, nz);
        T eh2 = ny * enz;
        T eh3 = nz * eny;
        T eh = eh1 + eh2 + eh3;

        std::array<T, 6> vec1y = {H, eh, -nx, -enx, -nx, -enx};
        std::array<T, 6> vec2y = {z0, z0, s, s, es, es};

        T y = dot_fma(vec1y, vec2y);

        T res_x = -x/S_2;
        T res_y = -y/S_2;
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
    static std::tuple<T, T> AccuDOP(T a, T b, T c, T d) {
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
    static std::tuple<T, T> CompDotC(const std::array<T, N>& v1, const std::array<T, N>& v2) {
        auto [s, c] = two_prod_fma(v1[0], v2[0]);
        for (std::size_t i = 1; i < N; ++i) {
            auto [p, pi] = two_prod_fma(v1[i], v2[i]);
            T sigma;
            std::tie(s, sigma) = two_sum(s, p); // Potential opportunity for speedup: do we really need this accuracy in the error term?
            c += pi + sigma;
        }
        return std::make_tuple(s, c);  
    }

    template <typename T, std::size_t N>
    static T dot_fma(const std::array<T, N>& v1, const std::array<T, N>& v2) {
        auto [s, c] = CompDotC(v1, v2);
        return s + c;
    }

    template <typename T>
    static std::tuple<T, T> sum_non_neg(T A_high, T A_low, T B_high, T B_low) {
        // High part of the sum using two_sum
        auto [H, h] = two_sum(A_high, B_high);

        // Low part of the sum
        T c = A_low + B_low;
        T d = h + c;


        return fast_two_sum(H, d);
    }

    template <typename T, std::size_t N>
    static std::tuple<T, T> SumOfSquaresC(const std::array<T, N>& vals, const std::array<T, N>& errs) {

        T S;
        T s;
        S = 0.0;
        s = 0.0;

        for (std::size_t i = 0; i < N; ++i) {
            auto [P, p] = two_prod_fma(vals[i], vals[i]); 
            std::tie(S, s) = sum_non_neg(S, s, P, p);  
        }


        T R = dot_fma(vals, errs);
        auto [S_new, s_err_new] = fast_two_sum(S, T(2 * R + s));

        return std::make_tuple(S_new, s_err_new); 
    }
};
#endif //CODES_GCACONSTLATINTERSECTIONS_H
