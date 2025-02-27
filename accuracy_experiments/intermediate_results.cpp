// intermediate_results.cpp
#include "Benchmarks.h"
#include "GCAConstLatIntersections.h"
#include "arcs_helpers.h"
#include "csv_helpers.h"
#include <iostream>
#include <boost/multiprecision/gmp.hpp>
#include <iostream>
#include <mp++/mp++.hpp>


using namespace boost::multiprecision;
template<typename T>
using V3_T = Eigen::Matrix<T, 3, 1>;
namespace bg = boost::geometry;
template<typename T>
using Arc_T = Eigen::Matrix<T, 3, 2>; // Matrix of 3 rows and 2 columns
namespace fs = std::filesystem;

void run_intermediate_result_exponent(int mpfr_precision,Arc_T<double> arc, double constZ,
                                  intermediate_values_our& iv_our, intermediate_values_compare<double>& ivc_double,
                                  intermediate_values_compare<mppp::real128>& ivc_quadruple,
                                  intermediate_values_compare<mpfr_float>& ivc_MPFR){
    AccurateIntersection ai;

    // Reconstruction of the pointA and pointB
    V3_T<double> pointA = arc.col(0);
    V3_T<double> pointB = arc.col(1);
    // convert pointA and pointB to array
    std::array<double,3> pointA_double = {pointA.x(), pointA.y(), pointA.z()};
    std::array<double,3> pointB_double = {pointB.x(), pointB.y(), pointB.z()};

    // Now calculate our final result
    std::tuple<double, double> intersection_ours = ai.gca_constLat_intersection_accurate(pointA_double, pointB_double, constZ,iv_our);


    // Now calculate the intersection using the double
    std::tuple<double,double> intersection_double = gca_constLat_intersection_coordinates_newEqn<double>(pointA, pointB, constZ, ivc_double);

    //Now calculate the intersection using the quadruple
    Eigen::Matrix<mppp::real128, 3, 1> pointA_real128 = pointA.cast<mppp::real128>();
    Eigen::Matrix<mppp::real128, 3, 1> pointB_real128 = pointB.cast<mppp::real128>();
    mppp::real128 constZ_real128 = static_cast<mppp::real128>(constZ);

    std::tuple<mppp::real128, mppp::real128> intersection_quadruple = gca_constLat_intersection_coordinates_newEqn<mppp::real128>(pointA_real128, pointB_real128, constZ_real128, ivc_quadruple);


    mpfr_float::default_precision(mpfr_precision);

    // Cast the pointA and pointB to mpfr_float
    V3_T<mpfr_float> pointA_mpfr(pointA.x(), pointA.y(), pointA.z());
    V3_T<mpfr_float> pointB_mpfr(pointB.x(), pointB.y(), pointB.z());

    // Now calculate the intersection using the MPFR
    std::tuple<mpfr_float, mpfr_float> intersection_MPFR = gca_constLat_intersection_coordinates_newEqn<mpfr_float>(pointA_mpfr, pointB_mpfr, constZ, ivc_MPFR);

}

void write_intermediate_result_exponent(int mpfr_precision,
                                        std::vector<intermediate_values_our>& ivs_our,
                                        std::vector<intermediate_values_compare<double>> ivcs_double,
                                        std::vector<intermediate_values_compare<mppp::real128>> ivcs_quadruple,
                                        std::vector<intermediate_values_compare<mpfr_float>> ivcs_MPFR,
                                        const std::string& filedirect,
                                        const std::string& filePrefix) {
    
    // Check if the directory exists, and if not, create it
    if (!fs::exists(filedirect)) {
        fs::create_directories(filedirect);  // Create the directory if it doesn't exist
    }

    // Now write the results to the file
    writeIV_Our(ivs_our, filedirect+filePrefix+"_intermediate_results_our.csv");
    writeIVC_Double(ivcs_double, filedirect+filePrefix+"_intermediate_results_double.csv");
    writeIVC_Quadruple(ivcs_quadruple, filedirect+filePrefix+"_intermediate_results_quadruple.csv");
    std::string filename = filedirect+filePrefix+"_intermediate_result_MPFR_" + std::to_string(mpfr_precision) + ".csv";
    writeIVC_MPFR(ivcs_MPFR, filename, mpfr_precision);

}

void run_write_intermediate_results_exponent(int mpfr_precision, 
                                             const std::vector<Arc_T<double>>& arcs, 
                                             const std::vector<double>& constZs,
                                             std::vector<intermediate_values_our>& ivs_our,
                                             std::vector<intermediate_values_compare<double>>& ivcs_double,
                                             std::vector<intermediate_values_compare<mppp::real128>>& ivcs_quadruple,
                                             std::vector<intermediate_values_compare<mpfr_float>>& ivcs_MPFR,
                                             const std::string& basePath,  // Add basePath as parameter
                                             const std::string& filePrefix) {
    // Run intermediate result calculations for each arc
    for (size_t i = 0; i < arcs.size(); ++i) {
        run_intermediate_result_exponent(mpfr_precision, arcs[i], constZs[i], ivs_our[i], ivcs_double[i], ivcs_quadruple[i], ivcs_MPFR[i]);
    }

    // Write the intermediate results to the specified directory
    write_intermediate_result_exponent(mpfr_precision, ivs_our, ivcs_double, ivcs_quadruple, ivcs_MPFR, basePath + "/intermediate_results/", filePrefix);
}



int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <mpfr_precisions> <latitudes>\n";
        return 1;
    }

    // Parse the MPFR precisions and latitudes from the input arguments
    auto mpfr_precisions = parseMPFRPrecisions(argv[1]);
    auto latitudes = parseDoubles(argv[2]);

    // Set the base path to the project root directory (assumed one level up from the current path)
    fs::path basePath = fs::current_path().parent_path();

    // Iterate over the MPFR precisions and latitude ranges
    for (auto mpfr_precision : mpfr_precisions) {
        for (size_t i = 0; i < latitudes.size() - 1; ++i) {
            // Process latitude range
            double startLat = latitudes[i];
            double endLat = latitudes[i + 1];

            // Generate the file names based on latitude ranges, formatted without decimal points
            std::string startLatStr = formatLatitude(startLat);
            std::string endLatStr = formatLatitude(endLat);
            std::string lat_range = startLatStr + "_" + endLatStr;
            std::string arc_file_name = lat_range + "Arcs_Exponent.csv";

            // Path to the arc file (from /generated_arcs directory)
            std::string arc_path = (basePath / "generated_arcs" / arc_file_name).string();

            // Read the region arcs and constant Z values from the CSV file
            std::vector<double> constZs;
            std::vector<Arc_T<double>> region_arcs = ReadGreatCircleArcsFromCSVExponent(arc_path, constZs);

            // Prepare vectors for intermediate results
            std::vector<intermediate_values_our> ivs_our(region_arcs.size());
            std::vector<intermediate_values_compare<double>> ivcs_double(region_arcs.size());
            std::vector<intermediate_values_compare<mppp::real128>> ivcs_quadruple(region_arcs.size());
            std::vector<intermediate_values_compare<mpfr_float>> ivcs_MPFR(region_arcs.size());

            // Run and write the intermediate results
            run_write_intermediate_results_exponent(mpfr_precision, region_arcs, constZs, ivs_our, ivcs_double, ivcs_quadruple, ivcs_MPFR, basePath.string(), lat_range);
        }
    }

    return 0;
}
