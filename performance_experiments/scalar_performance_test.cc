#include <tuple>
#include <iostream>
#include <chrono>
#include "simd_fma.hh"
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include "generate_input.hh"
#include <Eigen/Dense>
#include <utility>
#include <cmath>
constexpr size_t numTests = 100;
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

static std::tuple<double, double> new_implementation(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const double &constZ) {
    // Calculate the cross product, which gives us the vector n
    using T = double;
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


    // Calculate p_x and p_y, which are the x and y coordinates of the intersection
    T p_x = - (constZ * nx * nz + s_tilde * ny)/nx_squared_plus_ny_squared;

    T p_y = - (constZ * ny * nz - s_tilde * nx)/nx_squared_plus_ny_squared; // Only consider the minus

    //Print out the calculated p_x
//    std::cout << "p_x: " << p_x << std::endl;

    return {p_x, p_y};
}

static std::tuple<double, double> old_implementation(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const double constZ) {
    // Calculate the cross product, which gives us the vector n and normalize it
    Eigen::Vector3d n = simd_cross(pointA, pointB);
    n.normalize();

    // Extract nx, ny, nz components from vector n
    double nx = n[0];
    double ny = n[1];
    double nz = n[2];


    double nx_squared = nx * nx;
    double ny_squared = ny * ny;
    double nz_squared = nz * nz;
    double nx_squared_plus_ny_squared = nx_squared + ny_squared;


    double s = sqrt(nx_squared_plus_ny_squared - constZ * constZ);

    // Calculate p_x and p_y, which are the x and y coordinates of the intersection
    double p_x = - (constZ * nx * nz + s * ny) / nx_squared_plus_ny_squared;
    double p_y = - (constZ * ny * nz - s * nx) / nx_squared_plus_ny_squared; // Only consider the minus


    return {p_x, p_y};
}

int main(int argc, const char *argv[]) {
    static constexpr size_t dataSize = 1000000;
    Eigen::MatrixXd ptsA, ptsB;
    Eigen::VectorXd lattitudes;
    generateInput(dataSize, ptsA, ptsB, lattitudes);
    Eigen::VectorXd X(dataSize), Y(dataSize);

    // Run the benchmark and save the average times
    for (size_t repeat = 0; repeat < 1; ++repeat) {
        std::cout << "\n========= Naive Floating Point Results =========\n" << std::endl;
        // To store the timings for each method
        double time_new = 0.0, time_old = 0.0;


        // Measure time for the "new" method
        {
            auto start_new = std::chrono::high_resolution_clock::now();
            for (size_t i = 0; i < numTests; ++i) {
                Eigen::Vector3d x1, x2;
                double z0;
                
                for (size_t i = 0; i < dataSize; ++i) {
                    for (size_t c = 0; c < 3; ++c) {
                        x1(c) = ptsA(i, c);  // Assign directly for scalar T
                        x2(c) = ptsB(i, c);
                    }
                    z0 = lattitudes(i);

                    // Call the provided Intersection Function
                    auto [px, py] = new_implementation(x1, x2, z0);

                    X[i] = px;
                    Y[i] = py;
                }
            }
            std::chrono::duration<double> elapsed_new = std::chrono::high_resolution_clock::now() - start_new;
            time_new = elapsed_new.count();
        }

        // Measure time for the "old" method
        {
            auto start_old = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < numTests; ++i) {
            Eigen::Vector3d x1, x2;
            double z0;
            
            for (size_t i = 0; i < dataSize; ++i) {
                for (size_t c = 0; c < 3; ++c) {
                    x1(c) = ptsA(i, c);  // Assign directly for scalar T
                    x2(c) = ptsB(i, c);
                }
                z0 = lattitudes(i);

                // Call the provided Intersection Function
                auto [px, py] = old_implementation(x1, x2, z0);

                X[i] = px;
                Y[i] = py;
            }
        }
            std::chrono::duration<double> elapsed_old = std::chrono::high_resolution_clock::now() - start_old;
            time_old = elapsed_old.count();
        }


        double avg_time_new = time_new / (numTests * dataSize);
        double avg_time_old = time_old / (numTests * dataSize);

        // Print out the FP results
        std::cout<< " \n Naive FP old equation: " << avg_time_old <<" seconds \n" << std::endl;
        std::cout<< " \n Naive FP new equation: " << avg_time_new <<" seconds \n" << std::endl;
    }



    return 0;
}