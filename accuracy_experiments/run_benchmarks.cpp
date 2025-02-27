#include "Benchmarks.h"
#include "GCAConstLatIntersections.h"
#include "csv_helpers.h"
#include <string>
namespace fs = std::filesystem;


using namespace boost::multiprecision;
template<typename T>
using V3_T = Eigen::Matrix<T, 3, 1>;
template<typename T>
using Arc_T = Eigen::Matrix<T, 3, 2>; // Matrix of 3 rows and 2 columns




int main(int argc, char* argv[]){
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <mpfr_precisions> <latitudes> <offsets>\n";
        return 1;
    }

    auto mpfr_precisions = parseInts(argv[1]);
    auto latitudes = parseDoubles(argv[2]);

    std::string basePath = fs::current_path().parent_path().string();


    for (size_t i = 0; i < latitudes.size() - 1; ++i) {
        double startLat = latitudes[i];
        double endLat = latitudes[i + 1];

        // Generate the file names based on latitude ranges, formatted without decimal points.
        std::string startLatStr = formatLatitude(startLat);
        std::string endLatStr = formatLatitude(endLat);
        std::string lat_range = startLatStr + "_" + endLatStr;
        std::string arc_file_name = lat_range + "Arcs_Exponent.csv";

        // Generate arcs for the current latitude range.
        std::string arc_path = basePath+"/generated_arcs/" + arc_file_name;
        std::vector<double> constZs;
        std::vector<Arc_T<double>> region_arcs = ReadGreatCircleArcsFromCSVExponent(arc_path,constZs);

        std::vector<Arc_T<double>> sanitized_arcs;
        std::vector<double> sanitized_constZs;
        sanitized_arcs = region_arcs;
        sanitized_constZs = constZs;
        std::string output_file_directory = basePath+"/benchmark_results/"+lat_range+"_";


        // std::vector<double> p_x_baseline(region_arcs.size(), 0.0);
        // std::vector<double> p_y_baseline(region_arcs.size(), 0.0);
        


        // Sanitize each arc and write the results to a file.
        for (int i = 0; i < region_arcs.size(); ++i) {
            Eigen::Matrix<double, 3, 2> arc = region_arcs[i];
            // Get the start and end points of the arc
            V3_T<double> startPoint = arc.col(0);
            V3_T<double> endPoint = arc.col(1);

            // Get the constant latitude, which is the average of the start and end points/s z values
            double constZ = constZs[i];


            AlgorithmsAccuracyBenchMarks accuracyBenchmark16(startPoint, endPoint, constZ, 16, p_x_baseline[i], p_y_baseline[i]);
            AlgorithmsAccuracyBenchMarks accuracyBenchmark17(startPoint, endPoint, constZ, 17, p_x_baseline[i], p_y_baseline[i]);

            if (!containsNaNOrInf(accuracyBenchmark16.get_floating_point_baselineEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_floating_point_newEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_floating_point_oldEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_our_accurate_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_mpfr_newEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_mpfr_oldEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_mpfr_baselineEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_our_accurate_result()) &&
                !containsNaNOrInf(accuracyBenchmark17.get_mpfr_newEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark17.get_mpfr_oldEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark17.get_mpfr_baselineEqn_result())
                
                
                ) {


                sanitized_arcs.push_back(region_arcs[i]);
                sanitized_constZs.push_back(constZs[i]);
            }
        }

        // // Now write the sanitized arcs to a file.
        WriteGreatCircleArcsToCSVExponent(arc_path,sanitized_arcs, sanitized_constZs);

        std::cout<<"Done writing back filtered arcs \n";

        Now run the benchmarks for the sanitized arcs.

        std::string ourResultsPath = output_file_directory + "our_results_double.csv";
        std::string float64ResultsPath = output_file_directory + "float64_results_double.csv";
        std::string quadrupleResultsPath = output_file_directory + "quadruple_results_double.csv";

        // Now open new files for writing
        std::ofstream csvOur(ourResultsPath, std::ios::out | std::ios::trunc);
        std::ofstream csvFloat64(float64ResultsPath, std::ios::out | std::ios::trunc);
        std::ofstream csvQuadruple(quadrupleResultsPath, std::ios::out | std::ios::trunc);
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


        for (int i = 0; i < sanitized_arcs.size(); ++i) {
            Eigen::Matrix<double, 3, 2> arc = sanitized_arcs[i];
            V3_T<double> Point_A = arc.col(0);
            V3_T<double> Point_B = arc.col(1);
            double z_0 = sanitized_constZs[i];

            //Run double benchmarks
            std::tuple<double, double> float_baseline= gca_constLat_intersection_coordinates_baselineEqn<double>(Point_A, Point_B, z_0);
            std::tuple<double, double> float_new= gca_constLat_intersection_coordinates_newEqn<double>(Point_A, Point_B, z_0);
            std::tuple<double, double> float_old= gca_constLat_intersection_coordinates_oldEqn<double>(Point_A, Point_B, z_0);

            //Run quadruple benchmarks
            V3_T<mppp::real128> pointA_quad;
            V3_T<mppp::real128> pointB_quad;
            for (int i = 0; i < 3; ++i) {
                pointA_quad(i) = mppp::real128(Point_A(i));
                pointB_quad(i) = mppp::real128(Point_B(i));
            }
            mppp::real128 constZ_quad = mppp::real128(z_0);

            std::tuple<mppp::real128, mppp::real128> quadruple_baseline= gca_constLat_intersection_coordinates_baselineEqn<mppp::real128>(pointA_quad, pointB_quad, constZ_quad);
            std::tuple<mppp::real128, mppp::real128> quadruple_new= gca_constLat_intersection_coordinates_newEqn<mppp::real128>(pointA_quad, pointB_quad, constZ_quad);
            std::tuple<mppp::real128, mppp::real128> quadruple_old= gca_constLat_intersection_coordinates_oldEqn<mppp::real128>(pointA_quad, pointB_quad, constZ_quad);

            //Run our benchmarks
            std::array<double, 3> pointA_double = {Point_A(0), Point_A(1), Point_A(2)};
            std::array<double, 3> pointB_double = {Point_B(0), Point_B(1), Point_B(2)};
            AccurateIntersection gcaHelper;
            std::tuple<double, double> our_accurate= gcaHelper.gca_constLat_intersection_accurate(pointA_double, pointB_double, z_0);

            // Decompose and write float64 results
            try {
                DecomposedFloat xNew = createDecomposedFloat(std::get<0>(float_new), "get_floating_point_newEqn_result - x");
                DecomposedFloat yNew = createDecomposedFloat(std::get<1>(float_new), "get_floating_point_newEqn_result - y");
                DecomposedFloat xOld = createDecomposedFloat(std::get<0>(float_old), "get_floating_point_oldEqn_result - x");
                DecomposedFloat yOld = createDecomposedFloat(std::get<1>(float_old), "get_floating_point_oldEqn_result - y");
                DecomposedFloat xBaseline = createDecomposedFloat(std::get<0>(float_baseline), "get_floating_point_baselineEqn_result - x");
                DecomposedFloat yBaseline = createDecomposedFloat(std::get<1>(float_baseline), "get_floating_point_baselineEqn_result - y");
                if (csvFloat64.is_open()) {
                    csvFloat64 << xNew.significand << "," << xNew.exponent << ","
                            << yNew.significand << "," << yNew.exponent << ","
                            << xOld.significand << "," << xOld.exponent << ","
                            << yOld.significand << "," << yOld.exponent << ","
                            << xBaseline.significand << "," << xBaseline.exponent << ","
                            << yBaseline.significand << "," << yBaseline.exponent << "\n";
                }
            } catch (...) {}
            
            
            // Decompose and write our results
            try {
                DecomposedFloat xOur = createDecomposedFloat(std::get<0>(our_accurate), "get_our_accurate_result - x");
                DecomposedFloat yOur = createDecomposedFloat(std::get<1>(our_accurate), "get_our_accurate_result - y");
                if (csvOur.is_open()) {
                    csvOur << xOur.significand << "," << xOur.exponent << ","
                        << yOur.significand << "," << yOur.exponent << "\n";
                }
            } catch (...) {}



        try {

                double xNewQuadruple = static_cast<double>(std::get<0>(quadruple_new));
                double yNewQuadruple = static_cast<double>(std::get<1>(quadruple_new));
                double xOldQuadruple = static_cast<double>(std::get<0>(quadruple_old));
                double yOldQuadruple = static_cast<double>(std::get<1>(quadruple_old));
                double xBaselineQuadruple = static_cast<double>(std::get<0>(quadruple_baseline));
                double yBaselineQuadruple = static_cast<double>(std::get<1>(quadruple_baseline));

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

        }

        csvOur.close();
        csvFloat64.close();
        csvQuadruple.close();


        //Now run the MPFR benchmark
        for (auto mpfr_precision : mpfr_precisions) {
            std::string mpfrResultsPath = output_file_directory + "mpfr_" + std::to_string(mpfr_precision) + "_results.csv";
            std::ofstream csvMPFR(mpfrResultsPath, std::ios::out | std::ios::trunc);
            
            if (csvMPFR.is_open()) {
                csvMPFR << std::scientific << std::setprecision(mpfr_precision);
                csvMPFR << "x_newEqn_significand,x_newEqn_exponent,y_newEqn_significand,y_newEqn_exponent,x_oldEqn_significand,x_oldEqn_exponent,y_oldEqn_significand,y_oldEqn_exponent,x_baselineEqn_significand,x__baselineEqn_exponent,y__baselineEqn_significand,y__baselineEqn_exponent\n";
            }

            for (int j = 0; j < sanitized_arcs.size(); ++j) {
                Eigen::Matrix<double, 3, 2> arc = sanitized_arcs[j];
                V3_T<double> Point_A = arc.col(0);
                V3_T<double> Point_B = arc.col(1);
                double z_0 = sanitized_constZs[j];

                mpfr_float::default_precision(mpfr_precision);
                V3_T<mpfr_float> pointA_mpfr = Point_A.cast<mpfr_float>();
                V3_T<mpfr_float> pointB_mpfr = Point_B.cast<mpfr_float>();
                mpfr_float constZ_mpfr = z_0;

                std::tuple<mpfr_float , mpfr_float> mpfr_baseline = gca_constLat_intersection_coordinates_baselineEqn<mpfr_float>(pointA_mpfr, pointB_mpfr, constZ_mpfr);
                std::tuple<mpfr_float , mpfr_float> mpfr_old = gca_constLat_intersection_coordinates_oldEqn<mpfr_float>(pointA_mpfr, pointB_mpfr, constZ_mpfr);                
                std::tuple<mpfr_float , mpfr_float> mpfr_new = gca_constLat_intersection_coordinates_newEqn<mpfr_float>(pointA_mpfr, pointB_mpfr, constZ_mpfr);                


                try {

                    DecomposedFloat xNewMPFR = createDecomposedFloat(static_cast<double>(std::get<0>(mpfr_new)), "get_mpfr_newEqn_result - x");
                    DecomposedFloat yNewMPFR = createDecomposedFloat(static_cast<double>(std::get<1>(mpfr_new)), "get_mpfr_newEqn_result - y");
                    DecomposedFloat xOldMPFR = createDecomposedFloat(static_cast<double>(std::get<0>(mpfr_old)), "get_mpfr_oldEqn_result - x");
                    DecomposedFloat yOldMPFR = createDecomposedFloat(static_cast<double>(std::get<1>(mpfr_old)), "get_mpfr_oldEqn_result - y");
                    DecomposedFloat xBaselineMPFR = createDecomposedFloat(static_cast<double>(std::get<0>(mpfr_baseline)), "get_mpfr_baselineEqn_result - x");
                    DecomposedFloat yBaselineMPFR = createDecomposedFloat(static_cast<double>(std::get<1>(mpfr_baseline)), "get_mpfr_baselineEqn_result - y");

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
                } catch (...) {
                    std::cout<<" MPFR "<< mpfr_precision << "Has error \n";
                }
            }
            csvMPFR.close();

        }


        

    }


    // Parse the offsets
    std::vector<double> offsets =  parseDoubles(argv[4]);




    for (size_t i = 0; i < offsets.size() - 1; ++i) {
        double startOffset =  offsets[i];
        double endOffset =  offsets[i + 1];
        std::string startOffsetStr = formatOffset(startOffset);
        std::string endOffsetStr = formatOffset(endOffset);
        std::string offset_range = startOffsetStr + "_" + endOffsetStr;
        std::string arc_file_name = offset_range + "Extreme_Arcs_Exponent.csv";

        // Generate arcs for the current latitude range.
        std::string arc_path = basePath+"/generated_arcs/" + arc_file_name;
        std::vector<double> constZs;
        std::vector<Arc_T<double>> region_arcs = ReadGreatCircleArcsFromCSVExponent(arc_path,constZs);

        std::vector<Arc_T<double>> sanitized_arcs;
        std::vector<double> sanitized_constZs;
        // sanitized_arcs = region_arcs;
        // sanitized_constZs = constZs;
        std::string output_file_directory = basePath+"/benchmark_results/"+offset_range+"_Extreme_";


        std::vector<double> p_x_baseline(region_arcs.size(), 0.0);
        std::vector<double> p_y_baseline(region_arcs.size(), 0.0);
        


        // Sanitize each arc and write the results to a file.
        for (int i = 0; i < region_arcs.size(); ++i) {
            Eigen::Matrix<double, 3, 2> arc = region_arcs[i];
            // Get the start and end points of the arc
            V3_T<double> startPoint = arc.col(0);
            V3_T<double> endPoint = arc.col(1);

            // Get the constant latitude, which is the average of the start and end points/s z values
            double constZ = constZs[i];


            AlgorithmsAccuracyBenchMarks accuracyBenchmark16(startPoint, endPoint, constZ, 16, p_x_baseline[i], p_y_baseline[i]);
            AlgorithmsAccuracyBenchMarks accuracyBenchmark17(startPoint, endPoint, constZ, 17, p_x_baseline[i], p_y_baseline[i]);
            AlgorithmsAccuracyBenchMarks accuracyBenchmark30(startPoint, endPoint, constZ, 30, p_x_baseline[i], p_y_baseline[i]);
            AlgorithmsAccuracyBenchMarks accuracyBenchmark32(startPoint, endPoint, constZ, 32, p_x_baseline[i], p_y_baseline[i]);

            if (!containsNaNOrInf(accuracyBenchmark16.get_floating_point_baselineEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_floating_point_newEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_floating_point_oldEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_our_accurate_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_mpfr_newEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_mpfr_oldEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_mpfr_baselineEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark16.get_our_accurate_result()) &&
                !containsNaNOrInf(accuracyBenchmark17.get_mpfr_newEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark17.get_mpfr_oldEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark17.get_mpfr_baselineEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark30.get_mpfr_newEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark30.get_mpfr_oldEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark30.get_mpfr_baselineEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark32.get_mpfr_newEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark32.get_mpfr_oldEqn_result()) &&
                !containsNaNOrInf(accuracyBenchmark32.get_mpfr_baselineEqn_result())
                
                
                ) {


                sanitized_arcs.push_back(region_arcs[i]);
                sanitized_constZs.push_back(constZs[i]);
            }
        }

        // // Now write the sanitized arcs to a file.
        WriteGreatCircleArcsToCSVExponent(arc_path,sanitized_arcs, sanitized_constZs);

        std::cout<<"Done writing back filtered arcs \n";

        // Now run the benchmarks for the sanitized arcs.

        std::string ourResultsPath = output_file_directory + "our_results_double.csv";
        std::string float64ResultsPath = output_file_directory + "float64_results_double.csv";
        std::string quadrupleResultsPath = output_file_directory + "quadruple_results_double.csv";

        // Now open new files for writing
        std::ofstream csvOur(ourResultsPath, std::ios::out | std::ios::trunc);
        std::ofstream csvFloat64(float64ResultsPath, std::ios::out | std::ios::trunc);
        std::ofstream csvQuadruple(quadrupleResultsPath, std::ios::out | std::ios::trunc);
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


        for (int i = 0; i < sanitized_arcs.size(); ++i) {
            Eigen::Matrix<double, 3, 2> arc = sanitized_arcs[i];
            V3_T<double> Point_A = arc.col(0);
            V3_T<double> Point_B = arc.col(1);
            double z_0 = sanitized_constZs[i];

            //Run double benchmarks
            std::tuple<double, double> float_baseline= gca_constLat_intersection_coordinates_baselineEqn<double>(Point_A, Point_B, z_0);
            std::tuple<double, double> float_new= gca_constLat_intersection_coordinates_newEqn<double>(Point_A, Point_B, z_0);
            std::tuple<double, double> float_old= gca_constLat_intersection_coordinates_oldEqn<double>(Point_A, Point_B, z_0);

            //Run quadruple benchmarks
            V3_T<mppp::real128> pointA_quad;
            V3_T<mppp::real128> pointB_quad;
            for (int i = 0; i < 3; ++i) {
                pointA_quad(i) = mppp::real128(Point_A(i));
                pointB_quad(i) = mppp::real128(Point_B(i));
            }
            mppp::real128 constZ_quad = mppp::real128(z_0);

            std::tuple<mppp::real128, mppp::real128> quadruple_baseline= gca_constLat_intersection_coordinates_baselineEqn<mppp::real128>(pointA_quad, pointB_quad, constZ_quad);
            std::tuple<mppp::real128, mppp::real128> quadruple_new= gca_constLat_intersection_coordinates_newEqn<mppp::real128>(pointA_quad, pointB_quad, constZ_quad);
            std::tuple<mppp::real128, mppp::real128> quadruple_old= gca_constLat_intersection_coordinates_oldEqn<mppp::real128>(pointA_quad, pointB_quad, constZ_quad);

            //Run our benchmarks
            std::array<double, 3> pointA_double = {Point_A(0), Point_A(1), Point_A(2)};
            std::array<double, 3> pointB_double = {Point_B(0), Point_B(1), Point_B(2)};
            AccurateIntersection gcaHelper;
            std::tuple<double, double> our_accurate= gcaHelper.gca_constLat_intersection_accurate(pointA_double, pointB_double, z_0);

            // Decompose and write float64 results
            try {
                DecomposedFloat xNew = createDecomposedFloat(std::get<0>(float_new), "get_floating_point_newEqn_result - x");
                DecomposedFloat yNew = createDecomposedFloat(std::get<1>(float_new), "get_floating_point_newEqn_result - y");
                DecomposedFloat xOld = createDecomposedFloat(std::get<0>(float_old), "get_floating_point_oldEqn_result - x");
                DecomposedFloat yOld = createDecomposedFloat(std::get<1>(float_old), "get_floating_point_oldEqn_result - y");
                DecomposedFloat xBaseline = createDecomposedFloat(std::get<0>(float_baseline), "get_floating_point_baselineEqn_result - x");
                DecomposedFloat yBaseline = createDecomposedFloat(std::get<1>(float_baseline), "get_floating_point_baselineEqn_result - y");
                if (csvFloat64.is_open()) {
                    csvFloat64 << xNew.significand << "," << xNew.exponent << ","
                            << yNew.significand << "," << yNew.exponent << ","
                            << xOld.significand << "," << xOld.exponent << ","
                            << yOld.significand << "," << yOld.exponent << ","
                            << xBaseline.significand << "," << xBaseline.exponent << ","
                            << yBaseline.significand << "," << yBaseline.exponent << "\n";
                }
            } catch (...) {}
            
            
            // Decompose and write our results
            try {
                DecomposedFloat xOur = createDecomposedFloat(std::get<0>(our_accurate), "get_our_accurate_result - x");
                DecomposedFloat yOur = createDecomposedFloat(std::get<1>(our_accurate), "get_our_accurate_result - y");
                if (csvOur.is_open()) {
                    csvOur << xOur.significand << "," << xOur.exponent << ","
                        << yOur.significand << "," << yOur.exponent << "\n";
                }
            } catch (...) {}



        try {

                double xNewQuadruple = static_cast<double>(std::get<0>(quadruple_new));
                double yNewQuadruple = static_cast<double>(std::get<1>(quadruple_new));
                double xOldQuadruple = static_cast<double>(std::get<0>(quadruple_old));
                double yOldQuadruple = static_cast<double>(std::get<1>(quadruple_old));
                double xBaselineQuadruple = static_cast<double>(std::get<0>(quadruple_baseline));
                double yBaselineQuadruple = static_cast<double>(std::get<1>(quadruple_baseline));

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

        }

        csvOur.close();
        csvFloat64.close();
        csvQuadruple.close();


        //Now run the MPFR benchmark
        for (auto mpfr_precision : mpfr_precisions) {
            std::string mpfrResultsPath = output_file_directory + "mpfr_" + std::to_string(mpfr_precision) + "_results.csv";
            std::ofstream csvMPFR(mpfrResultsPath, std::ios::out | std::ios::trunc);
            
            if (csvMPFR.is_open()) {
                csvMPFR << std::scientific << std::setprecision(mpfr_precision);
                csvMPFR << "x_newEqn_significand,x_newEqn_exponent,y_newEqn_significand,y_newEqn_exponent,x_oldEqn_significand,x_oldEqn_exponent,y_oldEqn_significand,y_oldEqn_exponent,x_baselineEqn_significand,x__baselineEqn_exponent,y__baselineEqn_significand,y__baselineEqn_exponent\n";
            }

            for (int j = 0; j < sanitized_arcs.size(); ++j) {
                Eigen::Matrix<double, 3, 2> arc = sanitized_arcs[j];
                V3_T<double> Point_A = arc.col(0);
                V3_T<double> Point_B = arc.col(1);
                double z_0 = sanitized_constZs[j];

                mpfr_float::default_precision(mpfr_precision);
                V3_T<mpfr_float> pointA_mpfr = Point_A.cast<mpfr_float>();
                V3_T<mpfr_float> pointB_mpfr = Point_B.cast<mpfr_float>();
                mpfr_float constZ_mpfr = z_0;

                std::tuple<mpfr_float , mpfr_float> mpfr_baseline = gca_constLat_intersection_coordinates_baselineEqn<mpfr_float>(pointA_mpfr, pointB_mpfr, constZ_mpfr);
                std::tuple<mpfr_float , mpfr_float> mpfr_old = gca_constLat_intersection_coordinates_oldEqn<mpfr_float>(pointA_mpfr, pointB_mpfr, constZ_mpfr);                
                std::tuple<mpfr_float , mpfr_float> mpfr_new = gca_constLat_intersection_coordinates_newEqn<mpfr_float>(pointA_mpfr, pointB_mpfr, constZ_mpfr);                


                try {

                    DecomposedFloat xNewMPFR = createDecomposedFloat(static_cast<double>(std::get<0>(mpfr_new)), "get_mpfr_newEqn_result - x");
                    DecomposedFloat yNewMPFR = createDecomposedFloat(static_cast<double>(std::get<1>(mpfr_new)), "get_mpfr_newEqn_result - y");
                    DecomposedFloat xOldMPFR = createDecomposedFloat(static_cast<double>(std::get<0>(mpfr_old)), "get_mpfr_oldEqn_result - x");
                    DecomposedFloat yOldMPFR = createDecomposedFloat(static_cast<double>(std::get<1>(mpfr_old)), "get_mpfr_oldEqn_result - y");
                    DecomposedFloat xBaselineMPFR = createDecomposedFloat(static_cast<double>(std::get<0>(mpfr_baseline)), "get_mpfr_baselineEqn_result - x");
                    DecomposedFloat yBaselineMPFR = createDecomposedFloat(static_cast<double>(std::get<1>(mpfr_baseline)), "get_mpfr_baselineEqn_result - y");

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
                } catch (...) {
                    std::cout<<" MPFR "<< mpfr_precision << "Has error \n";
                }
            }
            csvMPFR.close();

        }


        

    }
        
    
}