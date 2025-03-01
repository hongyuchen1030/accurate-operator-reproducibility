#ifndef FORWARD_DECLARATIONS_HH
#define FORWARD_DECLARATIONS_HH
#include <Eigen/Dense>

template<size_t VEC_WIDTH>
using Vec = Eigen::Array<double, VEC_WIDTH, 1>;

template <typename T, bool Parallel = false>
void compute_values_EFT(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y);

template <typename T, bool Parallel = false>
void compute_values(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y);

template <typename T, bool Parallel = false>
double run_benchmark_EFT(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y);

template <typename T, bool Parallel = false>
double run_benchmark(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y);

#endif /* end of include guard: FORWARD_DECLARATIONS_HH */
