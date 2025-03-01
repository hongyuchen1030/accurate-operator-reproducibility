#include <Eigen/Dense>
#include <utility>
#include <cmath>

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

std::tuple<double, double> new_implementation(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const double constZ) {
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
