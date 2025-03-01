#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include "generate_input.hh"

constexpr size_t numTests = 100;

template<typename T, size_t N>
using VN_T = Eigen::Matrix<T, N, 1>;
template<typename T>
using V3_T = VN_T<T, 3>;

std::tuple<double, double> old_implementation(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const double constZ);
std::tuple<double, double> new_implementation(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const double constZ);

int main(int argc, const char *argv[]) {
    static constexpr size_t dataSize = 1000000;
    Eigen::MatrixXd ptsA, ptsB;
    Eigen::VectorXd lattitudes;
    generateInput(dataSize, ptsA, ptsB, lattitudes);
    Eigen::VectorXd X(dataSize), Y(dataSize);

    // Measure time for the "new" method
    {
        auto start = std::chrono::high_resolution_clock::now();
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
        std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
        std::cout << "old implementation timing: " << elapsed.count() / (numTests * dataSize) << std::endl;
    }

    {
        auto start = std::chrono::high_resolution_clock::now();
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
        std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
        std::cout << "new implementation timing: " << elapsed.count() / (numTests * dataSize) << std::endl;
    }
    
    return 0;
}

