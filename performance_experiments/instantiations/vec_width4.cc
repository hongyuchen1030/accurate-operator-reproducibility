#include "../Vectorized_Benchmark_template.hh"

template void  compute_values_EFT<Vec<4>, false>(const Eigen::MatrixXd &, const Eigen::MatrixXd &, const Eigen::VectorXd &, Eigen::VectorXd &, Eigen::VectorXd &);
template void      compute_values<Vec<4>, false>(const Eigen::MatrixXd &, const Eigen::MatrixXd &, const Eigen::VectorXd &, Eigen::VectorXd &, Eigen::VectorXd &);
template double run_benchmark_EFT<Vec<4>, false>(const Eigen::MatrixXd &, const Eigen::MatrixXd &, const Eigen::VectorXd &, Eigen::VectorXd &, Eigen::VectorXd &);
template double     run_benchmark<Vec<4>, false>(const Eigen::MatrixXd &, const Eigen::MatrixXd &, const Eigen::VectorXd &, Eigen::VectorXd &, Eigen::VectorXd &);
