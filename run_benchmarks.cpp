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


static void prepareBenchmarksToMathematica(const std::vector<Arc_T<double>> arcs,
                                           const std::vector<double> const_zs,
                                           const std::string &arcs_filename,
                                           const std::string &benchmark_filedirectory,
                                           const int mpfr_precision) {

    //const z is the average of the z values of the two endpoints of the arc
    std::vector<double> sanitized_constZs;
    // now run the benchmark

    auto sanitized_arcs = run_and_write_benchmark_results(arcs, const_zs, benchmark_filedirectory,sanitized_constZs, mpfr_precision);

    WriteGreatCircleArcsToCSVExponent(arcs_filename,sanitized_arcs, sanitized_constZs);


}

int main(int argc, char* argv[]){
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <mpfr_precisions> <latitudes>\n";
        return 1;
    }

    // Parse the MPFR precisions and latitudes from the input arguments
    auto mpfr_precisions = parseMPFRPrecisions(argv[1]);
    auto latitudes = parseLatitudes(argv[2]);

    // Set the base path to the project root directory (assumed one level up from the current path)
    fs::path basePath = fs::current_path().parent_path();

    // Ensure the benchmark_results directory exists
    fs::path benchmark_results_dir = basePath / "benchmark_results";
    if (!fs::exists(benchmark_results_dir)) {
        fs::create_directory(benchmark_results_dir);
    }

    // Iterate over the MPFR precisions and latitude ranges
    for (auto mpfr_precision : mpfr_precisions) {
        for (size_t i = 0; i < latitudes.size() - 1; ++i) {
            double startLat = latitudes[i];
            double endLat = latitudes[i + 1];

            // Generate the file names based on latitude ranges, formatted without decimal points
            std::string startLatStr = formatLatitude(startLat);
            std::string endLatStr = formatLatitude(endLat);
            std::string lat_range = startLatStr + "_" + endLatStr;
            std::string arc_file_name = lat_range + "Arcs_Exponent.csv";

            // Path to the arc file (from /generated_arcs directory)
            fs::path arc_path = basePath / "generated_arcs" / arc_file_name;

            // Read the region arcs and constant Z values from the CSV file
            std::vector<double> constZs;
            std::vector<Arc_T<double>> region_arcs = ReadGreatCircleArcsFromCSVExponent(arc_path.string(), constZs);

            // Prepare the path for the benchmark results
            std::string benchmark_path = (benchmark_results_dir / (lat_range + "_")).string();

            // Prepare data for Mathematica benchmarks
            prepareBenchmarksToMathematica(region_arcs, constZs, arc_path.string(), benchmark_path, mpfr_precision);
        }
    }

    return 0;
}