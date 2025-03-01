#ifndef GENERATE_INPUT_HH
#define GENERATE_INPUT_HH

#include <Eigen/Dense>

inline void generateInput(size_t dataSize, Eigen::MatrixXd &ptsA, Eigen::MatrixXd &ptsB, Eigen::VectorXd &z) {
    ptsA.resize(dataSize, 3);
    ptsB.resize(dataSize, 3);
       z.resize(dataSize);

    // Define some test points and constant Z using the provided values
    double pointA_x_significand = 5028374390644146;
    int pointA_x_exponent = -53;
    double pointA_y_significand = -7472957205960756;
    int pointA_y_exponent = -53;
    double pointA_z_significand = 6254432438282003;
    int pointA_z_exponent = -79;

    double pointB_x_significand = 5167685454902838;
    int pointB_x_exponent = -53;
    double pointB_y_significand = -7377307466399399;
    int pointB_y_exponent = -53;
    double pointB_z_significand = 4525606513452550;
    int pointB_z_exponent = -78;

    double constZ_significand = 7998403280412384;
    int constZ_exponent = -79;

    // Constructing PointA, PointB, and constZ
    Eigen::Vector3d pointA(
        pointA_x_significand * pow(2.0, pointA_x_exponent),
        pointA_y_significand * pow(2.0, pointA_y_exponent),
        pointA_z_significand * pow(2.0, pointA_z_exponent)
    );

    Eigen::Vector3d pointB(
        pointB_x_significand * pow(2.0, pointB_x_exponent),
        pointB_y_significand * pow(2.0, pointB_y_exponent),
        pointB_z_significand * pow(2.0, pointB_z_exponent)
    );

    double constZ = constZ_significand * pow(2.0, constZ_exponent);

    for (size_t i = 0; i < dataSize; ++i) {
        ptsA.row(i) = pointA;
        ptsB.row(i) = pointB;
    }

    z.setConstant(constZ);

    // Perturb all inputs but the last (so that the first result is for the exact input data above.
    ptsA.bottomRows(dataSize - 1) += 1e-12 * Eigen::MatrixXd::Random(dataSize - 1, ptsA.cols());
    ptsB.bottomRows(dataSize - 1) += 1e-12 * Eigen::MatrixXd::Random(dataSize - 1, ptsA.cols());
    z.tail(dataSize - 1) += 1e-12 * Eigen::VectorXd::Random(dataSize - 1);
}

#endif /* end of include guard: GENERATE_INPUT_HH */
