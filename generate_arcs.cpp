#include "arcs_helpers.h"
#include "csv_helpers.h"
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;



int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <number_of_arcs> <latitudes>\n";
        return 1;
    }
    int num_arcs = std::stoi(argv[1]);

    // // Parse latitudes
    // std::vector<double> latitudes = parseLatitudes(argv[2]);

    // // Loop through latitudes and generate arcs
    // for (size_t i = 0; i < latitudes.size() - 1; ++i) {
    //     double startLat = latitudes[i];
    //     double endLat = latitudes[i + 1];
    //     std::string startLatStr = formatLatitude(startLat);
    //     std::string endLatStr = formatLatitude(endLat);
    //     std::string lat_range = startLatStr + "_" + endLatStr;
    //     std::string arc_file_name = lat_range + "Arcs_Exponent.csv";

    //     // Define the relative path to the ./generated_arcs/ directory
    //     fs::path project_root = fs::current_path().parent_path(); // Get current working directory
    //     fs::path arc_directory = project_root / "generated_arcs"; // Relative path to generated_arcs/
    //     fs::create_directories(arc_directory); // Ensure the directory exists

    //     fs::path arc_path = arc_directory / arc_file_name; // Full path to the arc file

    //     std::vector<double> constZs;
    //     std::vector<Arc_T<double>> region_arcs = generateGreatCircleArcsMPFR(num_arcs, startLat, endLat, 0.0, 360.0, constZs);


    //     // Write the arcs to the CSV file
    //     WriteGreatCircleArcsToCSVExponent(arc_path.string(), region_arcs, constZs); // Use arc_path.string() for the file path
    // }

    // Now set up the experiment for the arc up cases
    // Start and end values
    double start = 0.01;
    double end = 1e-15;

    // Generate the offsets
    std::vector<double> offsets;
    for (double value = start; value >= end; value /= 10) {
        offsets.push_back(value);
    }

    // loop through offsets and generate arcs, we still using the format for the Lattitude cases
        for (size_t i = 0; i <  offsets.size() - 1; ++i) {
        double startOffset =  offsets[i];
        double endOffset =  offsets[i + 1];
        std::string startOffsetStr = formatOffset(startOffset);
        std::string endOffsetStr = formatOffset(endOffset);
        std::string offset_range = startOffsetStr + "_" + endOffsetStr;
        std::string arc_file_name = offset_range + "Extreme_Arcs_Exponent.csv";

        // Define the relative path to the ./generated_arcs/ directory
        fs::path project_root = fs::current_path().parent_path(); // Get current working directory
        fs::path arc_directory = project_root / "generated_arcs"; // Relative path to generated_arcs/
        fs::create_directories(arc_directory); // Ensure the directory exists

        fs::path arc_path = arc_directory / arc_file_name; // Full path to the arc file

        std::vector<double> constZs;
        std::vector<Arc_T<double>> region_arcs = generateArcsUpGreatCircleArcsMPFR(num_arcs, startOffset, endOffset, constZs);


        // Write the arcs to the CSV file
        WriteGreatCircleArcsToCSVExponent(arc_path.string(), region_arcs, constZs); // Use arc_path.string() for the file path
    }
    
    return 0;
}
