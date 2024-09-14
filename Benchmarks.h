//
// Created by hongyu chen on 4/8/24.
//

#ifndef CODES_BENCHMARKS_H
#define CODES_BENCHMARKS_H

#include <vector>
#include <cmath>
#include <iostream>
#include "GCAConstLatIntersections.h"
#include "csv_helpers.h"
#include <boost/multiprecision/mpfr.hpp>
#include <array>
#include <chrono>
#include <fstream>
#include <mp++/mp++.hpp>

// Include your templates here

using namespace boost::multiprecision;
namespace mp = boost::multiprecision;
using namespace std::chrono;

template<typename T, size_t N>
using VN_T = Eigen::Matrix<T, N, 1>;
template<typename T>
using V3_T = VN_T<T, 3>;
template<typename T>
using Arc_T = Eigen::Matrix<T, 3, 2>; // Matrix of 3 rows and 2 columns




class NewAlgorithmsTimeBenchMarks {

public:

    template<typename FloatType>
    void baseComparison(int iterations, V3_T<FloatType> pointA, V3_T<FloatType> pointB, FloatType constZ) {
        auto start = high_resolution_clock::now();

        double sum = 0;
        double x = 0;
        double y = 0;
        for (int i = 0; i < iterations; ++i) {
            std::tuple<FloatType, FloatType> res = gca_constLat_intersection_coordinates_newEqn<FloatType>(pointA, pointB, constZ);
            sum += static_cast<double>(std::get<0>(res)) + static_cast<double>(std::get<1>(res));
            x += static_cast<double>(std::get<0>(res));
            y += static_cast<double>(std::get<1>(res));
        }

        auto end = high_resolution_clock::now();
        duration<double, std::milli> elapsed = end - start;


        std::cout << "Total result: " << sum << std::endl;
        std::cout << "Total x: " << x << std::endl;
        std::cout << "Total y: " << y << std::endl;

        std::cout << "Average time for " << typeid(FloatType).name() << ": "
                  << elapsed.count() / iterations << " ms\n";
    }

    void benchmarkOurMethod(int iterations,  std::array<double, 3> pointA, std::array<double, 3> pointB, double z0) {
        AccurateIntersection gcaHelper;
        auto start = high_resolution_clock::now();
        double sum = 0;
        double x = 0;
        double y = 0;
        for (int i = 0; i < iterations; ++i) {
            std::tuple<double, double> res = gcaHelper.gca_constLat_intersection_accurate(pointA, pointB, z0);
            sum += std::get<0>(res) + std::get<1>(res);
            x += std::get<0>(res);
            y += std::get<1>(res);
        }


        auto end = high_resolution_clock::now();
        duration<double, std::milli> elapsed = end - start;

        std::cout << "Total result: " << sum << std::endl;
        std::cout << "Total x: " << x << std::endl;
        std::cout << "Total y: " << y << std::endl;

        std::cout << "Average time for our accurate method "
                  << (elapsed.count() / iterations)
                  << " ms" << std::endl;
    }

};

class OldAlgorithmsTimeBenchMarks {

public:

    template<typename FloatType>
    void baseComparison(int iterations, V3_T<FloatType> pointA, V3_T<FloatType> pointB, FloatType constZ) {
        auto start = high_resolution_clock::now();

        double sum = 0;
        double x = 0;
        double y = 0;
        for (int i = 0; i < iterations; ++i) {
            std::tuple<FloatType, FloatType> res = gca_constLat_intersection_coordinates_oldEqn<FloatType>(pointA, pointB, constZ);
            sum += static_cast<double>(std::get<0>(res)) + static_cast<double>(std::get<1>(res));
            x += static_cast<double>(std::get<0>(res));
            y += static_cast<double>(std::get<1>(res));
        }

        auto end = high_resolution_clock::now();
        duration<double, std::milli> elapsed = end - start;


        std::cout << "Total result: " << sum << std::endl;
        std::cout << "Total x: " << x << std::endl;
        std::cout << "Total y: " << y << std::endl;

        std::cout << "Average time for " << typeid(FloatType).name() << ": "
                  << elapsed.count() / iterations << " ms\n";
    }


};



struct BenchmarkResults {
    mpfr_float_1000 average_double_new_error;
    mpfr_float_1000 average_mpfr_new_error;
    mpfr_float_1000 average_our_error;
    mpfr_float_1000 average_double_old_error;
    mpfr_float_1000 average_mpfr_old_error;
};



class AlgorithmsAccuracyBenchMarks {

private:
    V3_T<double> Point_A;
    V3_T<double> Point_B;
    int mpfr_precision;
    double z_0;

public:
    // Constructor
    AlgorithmsAccuracyBenchMarks(const V3_T<double>& pointA,
                                    const V3_T<double>& pointB,
                                    const double z0,
                                    const int mpfr_prec, const double & px_baseline,
                                 const double& py_baseline)
            : Point_A(pointA), // Normalize pointA
              Point_B(pointB), // Normalize pointB
              z_0(z0), // Assign z0 to z_0
              mpfr_precision(mpfr_prec) // Assign mpfr_prec to mpfr_precision
    {
        // Print the baseline result
    }

    // If you want to provide getters for your private members:
    V3_T<double> getPointA() const { return Point_A; }
    V3_T<double> getPointB() const { return Point_B; }
    double getZ0() const { return z_0; }

    std::tuple<double, double> get_floating_point_baselineEqn_result() {
        // Note the use of 'z_0', which is the member variable, not 'z0'
        return gca_constLat_intersection_coordinates_baselineEqn<double>(Point_A, Point_B, z_0);
    }
    

    std::tuple<double, double> get_floating_point_newEqn_result() {
        // Note the use of 'z_0', which is the member variable, not 'z0'
        return gca_constLat_intersection_coordinates_newEqn<double>(Point_A, Point_B, z_0);
    }


    std::tuple<double, double> get_floating_point_oldEqn_result() {
        // Note the use of 'z_0', which is the member variable, not 'z0'
        return gca_constLat_intersection_coordinates_oldEqn<double>(Point_A, Point_B, z_0);
    }

    std::tuple<mpfr_float , mpfr_float> get_mpfr_newEqn_result() {
        mpfr_float::default_precision(mpfr_precision);
        V3_T<mpfr_float> pointA_mpfr = Point_A.cast<mpfr_float>();
        V3_T<mpfr_float> pointB_mpfr = Point_B.cast<mpfr_float>();
        mpfr_float constZ_mpfr = z_0;

        return gca_constLat_intersection_coordinates_newEqn<mpfr_float>(pointA_mpfr, pointB_mpfr, constZ_mpfr);
    }

    std::tuple<mpfr_float , mpfr_float> get_mpfr_oldEqn_result() {
        mpfr_float::default_precision(mpfr_precision);
        V3_T<mpfr_float> pointA_mpfr = Point_A.cast<mpfr_float>();
        V3_T<mpfr_float> pointB_mpfr = Point_B.cast<mpfr_float>();
        mpfr_float constZ_mpfr = z_0;

        return gca_constLat_intersection_coordinates_oldEqn<mpfr_float>(pointA_mpfr, pointB_mpfr, constZ_mpfr);
    }

    std::tuple<mpfr_float , mpfr_float> get_mpfr_baselineEqn_result() {
        mpfr_float::default_precision(mpfr_precision);
        V3_T<mpfr_float> pointA_mpfr = Point_A.cast<mpfr_float>();
        V3_T<mpfr_float> pointB_mpfr = Point_B.cast<mpfr_float>();
        mpfr_float constZ_mpfr = z_0;

        return gca_constLat_intersection_coordinates_baselineEqn<mpfr_float>(pointA_mpfr, pointB_mpfr, constZ_mpfr);
    }

    std::tuple<double, double> get_our_accurate_result(){
        AccurateIntersection gcaHelper;
        // Convert Point_A and Point_B to std::array<double, 3>
        std::array<double, 3> pointA_double = {Point_A(0), Point_A(1), Point_A(2)};
        std::array<double, 3> pointB_double = {Point_B(0), Point_B(1), Point_B(2)};
        return gcaHelper.gca_constLat_intersection_accurate(pointA_double, pointB_double, z_0);

    }

    std::tuple<mppp::real128, mppp::real128> get_quadruple_newEqn_result() {
        V3_T<mppp::real128> pointA_quad;
        V3_T<mppp::real128> pointB_quad;
        for (int i = 0; i < 3; ++i) {
            pointA_quad(i) = mppp::real128(Point_A(i));
            pointB_quad(i) = mppp::real128(Point_B(i));
        }
        mppp::real128 constZ_quad = mppp::real128(z_0);

        return gca_constLat_intersection_coordinates_newEqn<mppp::real128>(pointA_quad, pointB_quad, constZ_quad);
    }

    std::tuple<mppp::real128, mppp::real128> get_quadruple_oldEqn_result() {
        V3_T<mppp::real128> pointA_quad;
        V3_T<mppp::real128> pointB_quad;
        for (int i = 0; i < 3; ++i) {
            pointA_quad(i) = mppp::real128(Point_A(i));
            pointB_quad(i) = mppp::real128(Point_B(i));
        }
        mppp::real128 constZ_quad = mppp::real128(z_0);

        return gca_constLat_intersection_coordinates_oldEqn<mppp::real128>(pointA_quad, pointB_quad, constZ_quad);
    }

    std::tuple<mppp::real128, mppp::real128> get_quadruple_baselineEqn_result() {
        V3_T<mppp::real128> pointA_quad;
        V3_T<mppp::real128> pointB_quad;
        for (int i = 0; i < 3; ++i) {
            pointA_quad(i) = mppp::real128(Point_A(i));
            pointB_quad(i) = mppp::real128(Point_B(i));
        }
        mppp::real128 constZ_quad = mppp::real128(z_0);

        return gca_constLat_intersection_coordinates_baselineEqn<mppp::real128>(pointA_quad, pointB_quad, constZ_quad);
    }




};

// Generic template function to check for NaN and Inf values
template<typename T>
bool containsNaNOrInf(const std::tuple<T, T>& results) {
    return boost::math::isnan(std::get<0>(results)) || boost::math::isnan(std::get<1>(results)) ||
           boost::math::isinf(std::get<0>(results)) || boost::math::isinf(std::get<1>(results));
}

// Function to handle DecomposedFloat creation and error logging
template<typename T>
DecomposedFloat createDecomposedFloat(T v, const std::string& methodName) {
    try {
        return DecomposedFloat(v);
    } catch (const std::exception& e) {
        std::cerr << "Exception in method: " << methodName << ", Value: " << v << ", Error: " << e.what() << std::endl;
        throw; // Re-throw the exception to maintain behavior
    }
}


static std::vector<Eigen::Matrix<double, 3, 2>> run_and_write_benchmark_results(const std::vector<Eigen::Matrix<double, 3, 2>>& arcs,
                                               const std::vector<double> &const_Zs,
                                               const std::string &output_file_directory,
                                               std::vector<double>& sanitized_constZs,
                                               const int mpfr_precision){

    std::vector<Arc_T<double>> sanitized_arcs;

// Paths for the output files
    std::string ourResultsPath = output_file_directory + "our_results_double.csv";
    std::string float64ResultsPath = output_file_directory + "float64_results_double.csv";
    std::string quadrupleResultsPath = output_file_directory + "quadruple_results_double.csv";
    std::string mpfrResultsPath = output_file_directory + "mpfr_" + std::to_string(mpfr_precision) + "_results.csv";



// Now open new files for writing
    std::ofstream csvOur(ourResultsPath, std::ios::out | std::ios::trunc);
    std::ofstream csvFloat64(float64ResultsPath, std::ios::out | std::ios::trunc);
    std::ofstream csvQuadruple(quadrupleResultsPath, std::ios::out | std::ios::trunc);
    std::ofstream csvMPFR(mpfrResultsPath, std::ios::out | std::ios::trunc);

    // Write headers for the files
    if (csvOur.is_open()) {
        csvOur << "x_significand,x_exponent,y_significand,y_exponent\n";
    }
    if (csvFloat64.is_open()) {
        csvFloat64 << "x_newEqn_significand,x_newEqn_exponent,y_newEqn_significand,y_newEqn_exponent,x_oldEqn_significand,x_oldEqn_exponent,y_oldEqn_significand,y_oldEqn_exponent,x_baselineEqn_significand,x__baselineEqn_exponent,y__baselineEqn_significand,y__baselineEqn_exponent\n";
    }
    if (csvQuadruple.is_open()) {
        csvQuadruple << "x_newEqn_significand,x_newEqn_exponent,y_newEqn_significand,y_newEqn_exponent,x_oldEqn_significand,x_oldEqn_exponent,y_oldEqn_significand,y_oldEqn_exponent,x_baselineEqn_significand,x__baselineEqn_exponent,y__baselineEqn_significand,y__baselineEqn_exponent\n";
    }
    if (csvMPFR.is_open()) {
        csvMPFR << std::scientific << std::setprecision(mpfr_precision);
        csvMPFR << "x_newEqn_significand,x_newEqn_exponent,y_newEqn_significand,y_newEqn_exponent,x_oldEqn_significand,x_oldEqn_exponent,y_oldEqn_significand,y_oldEqn_exponent,x_baselineEqn_significand,x__baselineEqn_exponent,y__baselineEqn_significand,y__baselineEqn_exponent\n";
    }

    // generate a fake baseline intersection x andy such that it is all 0.0
    std::vector<double> p_x_baseline(arcs.size(), 0.0);
    std::vector<double> p_y_baseline(arcs.size(), 0.0);

    for (int i = 0; i < arcs.size(); ++i) {
        Eigen::Matrix<double, 3, 2> arc = arcs[i];
        // Get the start and end points of the arc
        V3_T<double> startPoint = arc.col(0);
        V3_T<double> endPoint = arc.col(1);

        // Get the constant latitude, which is the average of the start and end points/s z values
        double constZ = const_Zs[i];

        AlgorithmsAccuracyBenchMarks accuracyBenchmark(startPoint, endPoint, constZ, mpfr_precision, p_x_baseline[i], p_y_baseline[i]);
        if (!containsNaNOrInf(accuracyBenchmark.get_floating_point_baselineEqn_result()) &&
            !containsNaNOrInf(accuracyBenchmark.get_floating_point_newEqn_result()) &&
            !containsNaNOrInf(accuracyBenchmark.get_floating_point_oldEqn_result()) &&
            !containsNaNOrInf(accuracyBenchmark.get_our_accurate_result()) &&
            !containsNaNOrInf(accuracyBenchmark.get_mpfr_newEqn_result()) &&
            !containsNaNOrInf(accuracyBenchmark.get_mpfr_oldEqn_result()) &&
            !containsNaNOrInf(accuracyBenchmark.get_mpfr_baselineEqn_result())) {


            sanitized_arcs.push_back(arcs[i]);
            sanitized_constZs.push_back(const_Zs[i]);

            // Decompose and write our results
            try {
                DecomposedFloat xOur = createDecomposedFloat(std::get<0>(accuracyBenchmark.get_our_accurate_result()), "get_our_accurate_result - x");
                DecomposedFloat yOur = createDecomposedFloat(std::get<1>(accuracyBenchmark.get_our_accurate_result()), "get_our_accurate_result - y");
                if (csvOur.is_open()) {
                    csvOur << xOur.significand << "," << xOur.exponent << ","
                        << yOur.significand << "," << yOur.exponent << "\n";
                }
            } catch (...) {}

            // Decompose and write float64 results
            try {
                DecomposedFloat xNew = createDecomposedFloat(std::get<0>(accuracyBenchmark.get_floating_point_newEqn_result()), "get_floating_point_newEqn_result - x");
                DecomposedFloat yNew = createDecomposedFloat(std::get<1>(accuracyBenchmark.get_floating_point_newEqn_result()), "get_floating_point_newEqn_result - y");
                DecomposedFloat xOld = createDecomposedFloat(std::get<0>(accuracyBenchmark.get_floating_point_oldEqn_result()), "get_floating_point_oldEqn_result - x");
                DecomposedFloat yOld = createDecomposedFloat(std::get<1>(accuracyBenchmark.get_floating_point_oldEqn_result()), "get_floating_point_oldEqn_result - y");
                DecomposedFloat xBaseline = createDecomposedFloat(std::get<0>(accuracyBenchmark.get_floating_point_baselineEqn_result()), "get_floating_point_baselineEqn_result - x");
                DecomposedFloat yBaseline = createDecomposedFloat(std::get<1>(accuracyBenchmark.get_floating_point_baselineEqn_result()), "get_floating_point_baselineEqn_result - y");
                if (csvFloat64.is_open()) {
                    csvFloat64 << xNew.significand << "," << xNew.exponent << ","
                            << yNew.significand << "," << yNew.exponent << ","
                            << xOld.significand << "," << xOld.exponent << ","
                            << yOld.significand << "," << yOld.exponent << ","
                            << xBaseline.significand << "," << xBaseline.exponent << ","
                            << yBaseline.significand << "," << yBaseline.exponent << "\n";
                }
            } catch (...) {}

        try {
                auto newEqnResult = accuracyBenchmark.get_quadruple_newEqn_result();
                auto oldEqnResult = accuracyBenchmark.get_quadruple_oldEqn_result();
                auto baselineEqnResult = accuracyBenchmark.get_quadruple_baselineEqn_result();

                double xNewQuadruple = static_cast<double>(std::get<0>(newEqnResult));
                double yNewQuadruple = static_cast<double>(std::get<1>(newEqnResult));
                double xOldQuadruple = static_cast<double>(std::get<0>(oldEqnResult));
                double yOldQuadruple = static_cast<double>(std::get<1>(oldEqnResult));
                double xBaselineQuadruple = static_cast<double>(std::get<0>(baselineEqnResult));
                double yBaselineQuadruple = static_cast<double>(std::get<1>(baselineEqnResult));

                DecomposedFloat xNewQuadrupleDecomposed = createDecomposedFloat(xNewQuadruple, "get_quadruple_newEqn_result - x");
                DecomposedFloat yNewQuadrupleDecomposed = createDecomposedFloat(yNewQuadruple, "get_quadruple_newEqn_result - y");
                DecomposedFloat xOldQuadrupleDecomposed = createDecomposedFloat(xOldQuadruple, "get_quadruple_oldEqn_result - x");
                DecomposedFloat yOldQuadrupleDecomposed = createDecomposedFloat(yOldQuadruple, "get_quadruple_oldEqn_result - y");
                DecomposedFloat xBaselineQuadrupleDecomposed = createDecomposedFloat(xBaselineQuadruple, "get_quadruple_baselineEqn_result - x");
                DecomposedFloat yBaselineQuadrupleDecomposed = createDecomposedFloat(yBaselineQuadruple, "get_quadruple_baselineEqn_result - y");

                if (csvQuadruple.is_open()) {
                    csvQuadruple.precision(16);
                    csvQuadruple << std::scientific;
                    csvQuadruple << xNewQuadrupleDecomposed.significand << "," << xNewQuadrupleDecomposed.exponent << ","
                                << yNewQuadrupleDecomposed.significand << "," << yNewQuadrupleDecomposed.exponent << ","
                                << xOldQuadrupleDecomposed.significand << "," << xOldQuadrupleDecomposed.exponent << ","
                                << yOldQuadrupleDecomposed.significand << "," << yOldQuadrupleDecomposed.exponent << ","
                                << xBaselineQuadrupleDecomposed.significand << "," << xBaselineQuadrupleDecomposed.exponent << ","
                                << yBaselineQuadrupleDecomposed.significand << "," << yBaselineQuadrupleDecomposed.exponent << "\n";
                }
            } catch (...) {}


            // Decompose and write MPFR results
            try {
                DecomposedFloat xNewMPFR = createDecomposedFloat(static_cast<double>(std::get<0>(accuracyBenchmark.get_mpfr_newEqn_result())), "get_mpfr_newEqn_result - x");
                DecomposedFloat yNewMPFR = createDecomposedFloat(static_cast<double>(std::get<1>(accuracyBenchmark.get_mpfr_newEqn_result())), "get_mpfr_newEqn_result - y");
                DecomposedFloat xOldMPFR = createDecomposedFloat(static_cast<double>(std::get<0>(accuracyBenchmark.get_mpfr_oldEqn_result())), "get_mpfr_oldEqn_result - x");
                DecomposedFloat yOldMPFR = createDecomposedFloat(static_cast<double>(std::get<1>(accuracyBenchmark.get_mpfr_oldEqn_result())), "get_mpfr_oldEqn_result - y");
                DecomposedFloat xBaselineMPFR = createDecomposedFloat(static_cast<double>(std::get<0>(accuracyBenchmark.get_mpfr_baselineEqn_result())), "get_mpfr_baselineEqn_result - x");
                DecomposedFloat yBaselineMPFR = createDecomposedFloat(static_cast<double>(std::get<1>(accuracyBenchmark.get_mpfr_baselineEqn_result())), "get_mpfr_baselineEqn_result - y");

                if (csvMPFR.is_open()) {
                    csvMPFR.precision(32);
                    csvMPFR << std::scientific;
                    csvMPFR << xNewMPFR.significand << "," << xNewMPFR.exponent << ","
                            << yNewMPFR.significand << "," << yNewMPFR.exponent << ","
                            << xOldMPFR.significand << "," << xOldMPFR.exponent << ","
                            << yOldMPFR.significand << "," << yOldMPFR.exponent << ","
                            << xBaselineMPFR.significand << "," << xBaselineMPFR.exponent << ","
                            << yBaselineMPFR.significand << "," << yBaselineMPFR.exponent << "\n";
                }
            } catch (...) {}

        }
    }
    csvOur.close();
    csvFloat64.close();
    csvMPFR.close();

    return sanitized_arcs;

}




#endif //CODES_BENCHMARKS_H
