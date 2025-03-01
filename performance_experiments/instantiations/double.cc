#include "../Vectorized_Benchmark_template.hh"

template void  compute_values_EFT<double, false>(const Eigen::MatrixXd &, const Eigen::MatrixXd &, const Eigen::VectorXd &, Eigen::VectorXd &, Eigen::VectorXd &);
template void      compute_values<double, false>(const Eigen::MatrixXd &, const Eigen::MatrixXd &, const Eigen::VectorXd &, Eigen::VectorXd &, Eigen::VectorXd &);
template double run_benchmark_EFT<double, false>(const Eigen::MatrixXd &, const Eigen::MatrixXd &, const Eigen::VectorXd &, Eigen::VectorXd &, Eigen::VectorXd &);
template double     run_benchmark<double, false>(const Eigen::MatrixXd &, const Eigen::MatrixXd &, const Eigen::VectorXd &, Eigen::VectorXd &, Eigen::VectorXd &);
