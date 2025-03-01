#include <Eigen/Dense>
#include <utility>
#include <cmath>

template<typename T, size_t N>
using VN_T = Eigen::Matrix<T, N, 1>;
template<typename T>
using V3_T = VN_T<T, 3>;

template<typename T>
static V3_T<T> simd_cross(const V3_T<T>& v1, const V3_T<T>& v2) {
    V3_T<T> n;
    n(0) = v1.coeff(1) * v2.coeff(2) - v1.coeff(2) * v2.coeff(1);  // x = Ay * Bz - Az * By
    n(1) = v1.coeff(2) * v2.coeff(0) - v1.coeff(0) * v2.coeff(2);  // y = Az * Bx - Ax * Bz
    n(2) = v1.coeff(0) * v2.coeff(1) - v1.coeff(1) * v2.coeff(0);  // z = Ax * By - Ay * Bx
    return n;
}

std::tuple<double, double> old_implementation(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const double constZ) {
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
