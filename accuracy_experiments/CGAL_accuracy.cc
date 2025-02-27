
#include <string>
#include <Eigen/Dense>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/iterator.h>
namespace fs = std::filesystem;


using namespace boost::multiprecision;
template<typename T>
using V3_T = Eigen::Matrix<T, 3, 1>;
template<typename T>
using Arc_T = Eigen::Matrix<T, 3, 2>; // Matrix of 3 rows and 2 columns

using SK = CGAL::Exact_spherical_kernel_3;
using PT = typename SK::Point_3;

PT to_cgal(const Eigen::Ref<const Eigen::Vector3d> &p) { return PT(p[0], p[1], p[2]); }

// Convert a CGAL point type into Eigen (rounding each of its components
// to floating point). Note that in the case of an intersection point,
// the input components are algebraic numbers of degree 2.
template<typename CGALPt3>
Eigen::Vector3d to_eigen(const CGALPt3 &p) {
    Eigen::Vector3d result;

    result << CGAL::to_double(p.x()),
              CGAL::to_double(p.y()),
              CGAL::to_double(p.z());

    return result;
}


template<typename T>
struct PrecisionTraits;

template<>
struct PrecisionTraits<float> { static constexpr size_t p = 24; };  // 24 bits for the mantissa

template<>
struct PrecisionTraits<double> { static constexpr size_t p = 53; }; // 53 bits for the mantissa
// Function to log errors
void logError(const std::string& methodName, const std::string& message, double value, int64_t significand, int64_t exponent) {
    std::cerr << "Error in method: " << methodName << ", Value: " << value << ", Message: " << message 
              << ", Significand: " << significand << ", Exponent: " << exponent << std::endl;
}
struct DecomposedFloat {
    long significand;
    long exponent;

    DecomposedFloat() : significand(0), exponent(0) {}

    // Constructor to directly set significand and exponent
    DecomposedFloat(long sig, long exp) : significand(sig), exponent(exp) {}

    template<typename T>
    DecomposedFloat(T v) {
        if (v == 0) {
            significand = 0;
            exponent = 0;
        } else {
            exponent = std::ilogb(v) - (std::numeric_limits<T>::digits - 1);
            T value_exp = std::ldexp(v, -exponent);
            significand = long(value_exp);
            if (significand == 0) {
                exponent = 0;
            }
            if (T(significand) != value_exp) {
                logError("DecomposedFloat", "Integer part was not integer!", v, significand, exponent);
                throw std::runtime_error("Integer part was not integer!");
            }
            if (toFloat<T>() != v) {
                logError("DecomposedFloat", "Could not reconstruct the original value!", v, significand, exponent);
                throw std::runtime_error("Could not reconstruct the original value!");
            }
        }
    }

    // Output stream operator
    friend std::ostream& operator<<(std::ostream& os, const DecomposedFloat &df) {
        return os << df.significand << " * 2**" << df.exponent;
    }

    template<typename T>
    T toFloat() const {
        return std::ldexp(T(significand), exponent);
    }

    // Convenience function to directly output as CSV
    std::string toCSV() const {
        return std::to_string(significand) + "," + std::to_string(exponent);
    }
};

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


std::string formatLatitude(double latitude) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(6) << latitude;
    std::string latStr = stream.str();
    for (auto& c : latStr) {
        if (c == '.') {
            c = '_';
        }
    }
    return latStr;
}

std::string formatOffset(double offset) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(16) << offset;
    std::string latStr = stream.str();
    for (auto& c : latStr) {
        if (c == '.') {
            c = '_';
        }
    }
    return latStr;
}


std::vector<double> parseDoubles(const std::string& input) {
    std::vector<double> latitudes;
    std::stringstream ss(input);
    std::string item;
    while (getline(ss, item, ',')) {
        latitudes.push_back(std::atof(item.c_str()));
    }
    return latitudes;
}

// Helper function to reconstruct double from significand and exponent
double reconstructDouble(long significand, long exponent) {
    return std::ldexp(static_cast<double>(significand), exponent);
}

// Function to read great circle arcs from a CSV file into a vector
std::vector<Arc_T<double>> ReadGreatCircleArcsFromCSVExponent(const std::string& filename, std::vector<double>& constZs) {
    std::vector<Arc_T<double>> arcs;
    std::ifstream csvFile(filename);
    if (!csvFile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return arcs;
    }

    std::string line;
    std::getline(csvFile, line); // Skip the header line

    while (std::getline(csvFile, line)) {
        std::stringstream ss(line);
        std::string value;
        std::vector<long> significands;
        std::vector<int> exponents;

        // Extract each comma-separated value and separate significands and exponents
        while (std::getline(ss, value, ',')) {
            long significand;
            int exponent;
            static bool toggle = true; // Toggle between reading significand and exponent

            if (toggle) { // Read significand
                significand = std::stol(value);
                significands.push_back(significand);
            } else { // Read exponent
                exponent = std::stoi(value);
                exponents.push_back(exponent);
            }
            toggle = !toggle; // Flip toggle for next read
        }


        // Reconstruct doubles and create an Arc_T
        Arc_T<double> arc;
        for (int i = 0; i < 6; ++i) {
            int row = i % 3;
            int col = i / 3;
            arc(row, col)= reconstructDouble(significands[i], exponents[i]);
        }
        arcs.push_back(arc);

        // Reconstruct constZ and add to the vector
        double constZ = reconstructDouble(significands[6], exponents[6]);
        constZs.push_back(constZ);
    }


    csvFile.close();
    return arcs;
}


int main(int argc, char* argv[]){
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <latitudes>\n";
        return 1;
    }

    auto latitudes = parseDoubles(argv[1]);

    std::string basePath = fs::current_path().parent_path().string();


    for (size_t i = 0; i < latitudes.size() - 1; ++i) {
        double startLat = latitudes[i];
        double endLat = latitudes[i + 1];

        // Generate the file names based on latitude ranges, formatted without decimal points.
        std::string startLatStr = formatLatitude(startLat);
        std::string endLatStr = formatLatitude(endLat);
        std::string lat_range = startLatStr + "_" + endLatStr;
        std::string arc_file_name = lat_range + "Arcs_Exponent.csv";

        // The arcs are santitized already
        std::string arc_path = basePath+"/generated_arcs/" + arc_file_name;
        std::vector<double> sanitized_constZs;
        std::vector<Arc_T<double>> sanitized_arcs = ReadGreatCircleArcsFromCSVExponent(arc_path,sanitized_constZs);

        
        std::string output_file_directory = basePath+"/benchmark_results/"+lat_range+"_";
        std::string cgalResultsPath = output_file_directory + "cgal_double.csv";
        std::ofstream csvCgal(cgalResultsPath, std::ios::out | std::ios::trunc);
        if (csvCgal.is_open()) {
            csvCgal << "x1_significand,x1_exponent,y1_significand,y1_exponent,x2_significand,x2_exponent,y2_significand,y2_exponent\n";
        }


        PT origin(0, 0, 0);
        SK::Vector_3 z_axis(0, 0, 1);
        SK::Sphere_3 unit_sphere(origin, 1);


        // unsigned integer indicates multiplicity of intersection point
        using Point_and_multiplicity = std::pair<SK::Circular_arc_point_3, unsigned>;

        // Now run the benchmarks for the sanitized arcs.
        for (int i = 0; i < sanitized_arcs.size(); ++i) {
            Eigen::Matrix<double, 3, 2> arc = sanitized_arcs[i];
            V3_T<double> Point_A = arc.col(0);
            V3_T<double> Point_B = arc.col(1);
            double z_0 = sanitized_constZs[i];

            // Great circle passing through A and B.
            SK::Circle_3 c1(unit_sphere, SK::Plane_3(origin, to_cgal(Point_A), to_cgal(Point_B)));

            // Line of constant lattiude
            SK::Circle_3 c2(unit_sphere, SK::Plane_3(PT(0, 0, z_0), z_axis));

            SK::Intersect_3 inter;
            std::vector<Point_and_multiplicity> intersections;

            // Only recover point-valued intersections (not circles)
            inter(c1, c2, CGAL::dispatch_or_drop_output<Point_and_multiplicity>(std::back_inserter(intersections)));
            
            double x_1 = CGAL::to_double(intersections[0].first.x());
            double y_1 = CGAL::to_double(intersections[0].first.y());

            double x_2 = CGAL::to_double(intersections[1].first.x());
            double y_2 = CGAL::to_double(intersections[1].first.y());

            // Now write the two intersection point into the csv file
                        // Decompose and write float64 results
            try {
                DecomposedFloat x1_double = createDecomposedFloat(x_1, "CGAL x1");
                DecomposedFloat y1_double = createDecomposedFloat(y_1, "CGAL y1");
                DecomposedFloat x2_double = createDecomposedFloat(x_2, "CGAL x2");
                DecomposedFloat y2_double = createDecomposedFloat(y_2, "CGAL y2");
  
                if (csvCgal.is_open()) {
                    csvCgal << x1_double.significand << "," << x1_double.exponent << ","
                            << y1_double.significand << "," << y1_double.exponent << ","
                            << x2_double.significand << "," << x2_double.exponent << ","
                            << y2_double.significand << "," << y2_double.exponent << "\n";
                }
            } catch (...) {}
            


        }
        csvCgal.close();




        

    }

    // Now run the extreme arcup case
    std::vector<double> offsets =  parseDoubles(argv[2]);
    for (size_t i = 0; i < offsets.size() - 1; ++i) {

        double startOffset =  offsets[i];
        double endOffset =  offsets[i + 1];
        std::string startOffsetStr = formatOffset(startOffset);
        std::string endOffsetStr = formatOffset(endOffset);
        std::string offset_range = startOffsetStr + "_" + endOffsetStr;
        std::string arc_file_name = offset_range + "Extreme_Arcs_Exponent.csv";

        // Generate arcs for the current latitude range.
        std::string arc_path = basePath+"/generated_arcs/" + arc_file_name;
        std::vector<double> sanitized_constZs;
        std::vector<Arc_T<double>> sanitized_arcs = ReadGreatCircleArcsFromCSVExponent(arc_path,sanitized_constZs);

        
        std::string output_file_directory = basePath+"/benchmark_results/"+offset_range+"_Extreme_";
        std::string cgalResultsPath = output_file_directory + "cgal_double.csv";
        std::ofstream csvCgal(cgalResultsPath, std::ios::out | std::ios::trunc);
        if (csvCgal.is_open()) {
            csvCgal << "x1_significand,x1_exponent,y1_significand,y1_exponent,x2_significand,x2_exponent,y2_significand,y2_exponent\n";
        }


        PT origin(0, 0, 0);
        SK::Vector_3 z_axis(0, 0, 1);
        SK::Sphere_3 unit_sphere(origin, 1);


        // unsigned integer indicates multiplicity of intersection point
        using Point_and_multiplicity = std::pair<SK::Circular_arc_point_3, unsigned>;

        // Now run the benchmarks for the sanitized arcs.
        for (int i = 0; i < sanitized_arcs.size(); ++i) {
            Eigen::Matrix<double, 3, 2> arc = sanitized_arcs[i];
            V3_T<double> Point_A = arc.col(0);
            V3_T<double> Point_B = arc.col(1);
            double z_0 = sanitized_constZs[i];

            // Great circle passing through A and B.
            SK::Circle_3 c1(unit_sphere, SK::Plane_3(origin, to_cgal(Point_A), to_cgal(Point_B)));

            // Line of constant lattiude
            SK::Circle_3 c2(unit_sphere, SK::Plane_3(PT(0, 0, z_0), z_axis));

            SK::Intersect_3 inter;
            std::vector<Point_and_multiplicity> intersections;

            // Only recover point-valued intersections (not circles)
            inter(c1, c2, CGAL::dispatch_or_drop_output<Point_and_multiplicity>(std::back_inserter(intersections)));
            
            double x_1 = CGAL::to_double(intersections[0].first.x());
            double y_1 = CGAL::to_double(intersections[0].first.y());

            double x_2 = CGAL::to_double(intersections[1].first.x());
            double y_2 = CGAL::to_double(intersections[1].first.y());

            // Now write the two intersection point into the csv file
                        // Decompose and write float64 results
            try {
                DecomposedFloat x1_double = createDecomposedFloat(x_1, "CGAL x1");
                DecomposedFloat y1_double = createDecomposedFloat(y_1, "CGAL y1");
                DecomposedFloat x2_double = createDecomposedFloat(x_2, "CGAL x2");
                DecomposedFloat y2_double = createDecomposedFloat(y_2, "CGAL y2");
  
                if (csvCgal.is_open()) {
                    csvCgal << x1_double.significand << "," << x1_double.exponent << ","
                            << y1_double.significand << "," << y1_double.exponent << ","
                            << x2_double.significand << "," << x2_double.exponent << ","
                            << y2_double.significand << "," << y2_double.exponent << "\n";
                }
            } catch (...) {}
            


        }
        csvCgal.close();




        

    }


}