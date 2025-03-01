#include <Eigen/Dense>
#include <utility>
#include <cmath>

template<typename T, size_t N>
using VN_T = Eigen::Matrix<T, N, 1>;
template<typename T>
using V3_T = VN_T<T, 3>;

// template<typename T>
// static V3_T<T> simd_cross(const V3_T<T>& v1, const V3_T<T>& v2) {
//     V3_T<T> n;
//     n(0) = v1.coeff(1) * v2.coeff(2) - v1.coeff(2) * v2.coeff(1);
//     n(1) = v1.coeff(2) * v2.coeff(0) - v1.coeff(0) * v2.coeff(2);
//     n(2) = v1.coeff(0) * v2.coeff(1) - v1.coeff(1) * v2.coeff(0);
//     return n;
// }

static std::tuple<double, double> new_implementation(const Eigen::Vector3d& pointA, const Eigen::Vector3d& pointB, const double constZ) {
    V3_T<double> n = simd_cross(pointA, pointB);
    double nx = n[0], ny = n[1], nz = n[2];
    double nx_squared_plus_ny_squared = nx * nx + ny * ny;
    double norm_n_squared = nx_squared_plus_ny_squared + nz * nz;
    double s_tilde = sqrt(nx_squared_plus_ny_squared - norm_n_squared * constZ * constZ);
    double p_x = -(constZ * nx * nz + s_tilde * ny) / nx_squared_plus_ny_squared;
    double p_y = -(constZ * ny * nz - s_tilde * nx) / nx_squared_plus_ny_squared;
    return {p_x, p_y};
}
