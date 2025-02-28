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


bool isValidArc(std::vector<AlgorithmsAccuracyBenchMarks>& benchmarks) {
    for (auto& benchmark : benchmarks) {
        if (containsNaNOrInf(benchmark.get_floating_point_baselineEqn_result()) ||
            containsNaNOrInf(benchmark.get_floating_point_newEqn_result()) ||
            containsNaNOrInf(benchmark.get_floating_point_oldEqn_result()) ||
            containsNaNOrInf(benchmark.get_our_accurate_result()) ||
            containsNaNOrInf(benchmark.get_mpfr_newEqn_result()) ||
            containsNaNOrInf(benchmark.get_mpfr_oldEqn_result()) ||
            containsNaNOrInf(benchmark.get_mpfr_baselineEqn_result())) {
            return false;
        }
    }
    return true;
}

void processArcs(const std::vector<double>& values, const std::string& arc_prefix, 
    const std::string& output_prefix, const std::vector<int>& mpfr_precisions, bool isLatitude) {

    std::string basePath = fs::current_path().parent_path().string();

    for (size_t i = 0; i < values.size() - 1; ++i) {
    double startVal = values[i];
    double endVal = values[i + 1];
    
    std::string startStr, endStr;
    if (isLatitude) {
        // Formatting for latitudes
        startStr = formatLatitude(startVal);
        endStr = formatLatitude(endVal);
    } else {
        // Formatting for offsets
        startStr = formatOffset(startVal);
        endStr = formatOffset(endVal);
    }
    std::string range = startStr + "_" + endStr;
    std::string arc_file_name = range + arc_prefix;

    std::string arc_path = basePath + "/generated_arcs/" + arc_file_name;
    std::vector<double> constZs;
    std::vector<Arc_T<double>> region_arcs = ReadGreatCircleArcsFromCSVExponent(arc_path, constZs);

    std::vector<Arc_T<double>> sanitized_arcs;
    std::vector<double> sanitized_constZs;
    std::string output_file_directory = basePath + "/benchmark_results/" + range + "_";

    std::vector<double> p_x_baseline(region_arcs.size(), 0.0);
    std::vector<double> p_y_baseline(region_arcs.size(), 0.0);

    // Process each arc
    for (size_t i = 0; i < region_arcs.size(); ++i) {
    Eigen::Matrix<double, 3, 2> arc = region_arcs[i];
    V3_T<double> startPoint = arc.col(0);
    V3_T<double> endPoint = arc.col(1);
    double constZ = constZs[i];

    std::vector<AlgorithmsAccuracyBenchMarks> benchmarks;
    for (int precision : mpfr_precisions) {
    benchmarks.emplace_back(startPoint, endPoint, constZ, precision, p_x_baseline[i], p_y_baseline[i]);
    }

    if (isValidArc(benchmarks)) {
    sanitized_arcs.push_back(region_arcs[i]);
    sanitized_constZs.push_back(constZs[i]);
    }
    }

    // Write the sanitized arcs to a file
    WriteGreatCircleArcsToCSVExponent(arc_path, sanitized_arcs, sanitized_constZs);
    std::cout << "Done writing back filtered arcs: " << arc_path << "\n";
    }
}


int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <mpfr_precisions> <latitudes> <offsets>\n";
        return 1;
    }

    // Dynamically parse MPFR precisions from input
    std::vector<int> mpfr_precisions = parseInts(argv[1]);
    std::vector<double> latitudes = parseDoubles(argv[2]);
    std::vector<double> offsets = parseDoubles(argv[3]);

    // Process latitude-based arcs
    processArcs(latitudes, "Arcs_Exponent.csv", "benchmark_results/", mpfr_precisions,true);

    // Process offset-based arcs
    processArcs(offsets, "ArcUp_Arcs_Exponent.csv", "benchmark_results/", mpfr_precisions,false);
}
