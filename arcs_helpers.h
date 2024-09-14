//
// Created by hongyu chen on 4/15/24.
//

#ifndef CODES_ARCS_GENERATOR_H
#define CODES_ARCS_GENERATOR_H
#include "Benchmarks.h"
#include "GCAConstLatIntersections.h"
#include <iostream>
#include <boost/multiprecision/gmp.hpp>
#include <iomanip> // for std::setprecision
#include <boost/geometry.hpp>
#include <Eigen/Dense>
#include <random>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/random.hpp>


using namespace boost::multiprecision;
template<typename T>
using V3_T = Eigen::Matrix<T, 3, 1>;
template<typename T>
using Arc_T = Eigen::Matrix<T, 3, 2>; // Matrix of 3 rows and 2 columns
namespace bg = boost::geometry;


// Function to convert degrees to radians
template<typename T>
static T degreesToRadians(T degrees) {
    const T pi = boost::math::constants::pi<T>();
    return degrees * (pi / 180.0);
}

// Templated function to convert latitude and longitude to Cartesian coordinates on a unit sphere
template<typename T>
static V3_T<T> latLonToXYZ(T latitude, T longitude) {
    T latRad = degreesToRadians(latitude);
    T lonRad = degreesToRadians(longitude);


    return V3_T<T>(
            cos(latRad) * cos(lonRad), // x
            cos(latRad) * sin(lonRad), // y
            sin(latRad)                // z
    ).normalized();
}

// Helper function to check the span against the desired span within a tolerance
template<typename T>
static bool isSpanCorrect(const V3_T<T>& v1, const V3_T<T>& v2, T desiredSpanRadians, double toleranceRadians) {
    T angleRadians = acos(v1.normalized().dot(v2.normalized()));
    return abs(angleRadians - desiredSpanRadians) <= toleranceRadians;
}


// Function to generate random great circle arcs with a specified span
static std::vector<Arc_T<double>> generateGreatCircleArcs(
        int numArcs,
        double lat_min, double lat_max,
        double lon_min, double lon_max,
        double spanDegrees = -1.0, // Optional: desired span in degrees
        double toleranceDegrees = 0.1)
{
    using namespace boost::random;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> latDist(lat_min, lat_max);
    std::uniform_real_distribution<double> lonDist(lon_min, lon_max);

    std::vector<Arc_T<double>> arcs;
    double spanRadians = degreesToRadians(spanDegrees);
    double toleranceRadians = degreesToRadians(toleranceDegrees);

    for (int i = 0; i < numArcs; ++i) {
        V3_T<double> startPoint, endPoint;
        do {
            double startLat = latDist(gen);
            double startLon = lonDist(gen);
            double endLat = latDist(gen);
            double endLon = lonDist(gen);

            startPoint = latLonToXYZ(startLat, startLon);
            endPoint = latLonToXYZ(endLat, endLon);
        } while (spanDegrees > 0 && !isSpanCorrect(startPoint, endPoint, spanRadians, toleranceRadians));

        Arc_T<double> arc;
        arc.col(0) = startPoint;
        arc.col(1) = endPoint;

        arcs.push_back(arc);
    }

    return arcs;
}


// Function to generate random great circle arcs with a specified span using MPFR precision
// The return will be casted to double precision
static std::vector<Arc_T<double>> generateGreatCircleArcsMPFR(
        int numArcs,
        mpfr_float_1000 lat_min, mpfr_float_1000 lat_max,
        mpfr_float_1000 lon_min, mpfr_float_1000 lon_max,
        std::vector<double>&constZs,
        mpfr_float_1000 spanDegrees = -1.0) // Optional: desired span in degrees)
{
//    std::random_device rd;
    std::mt19937 gen(0);
    std::uniform_real_distribution<double> Dist(0.0, 1.0);

    std::vector<Arc_T<double>> arcs;

    mpfr_float_1000 spanRadians = degreesToRadians(spanDegrees);


    for (int i = 0; i < numArcs; ++i) {
        V3_T<mpfr_float_1000> startPoint, endPoint;
        do {
            // Generate random points within 0.0 and 1.0
            double startLat_double = Dist(gen);
            double startLon_double = Dist(gen);
            double endLat_double = Dist(gen);
            double endLon_double = Dist(gen);

            // Scale the random points to the desired range
            mpfr_float_1000 startLat = lat_min + startLat_double * (lat_max - lat_min);
            mpfr_float_1000 startLon = lon_min + startLon_double * (lon_max - lon_min);
            mpfr_float_1000 endLat = lat_min + endLat_double * (lat_max - lat_min);
            mpfr_float_1000 endLon = lon_min + endLon_double * (lon_max - lon_min);

            startPoint = latLonToXYZ<mpfr_float_1000>(startLat, startLon);
            endPoint = latLonToXYZ<mpfr_float_1000>(endLat, endLon);
        } while (spanDegrees > 0 && !isSpanCorrect<mpfr_float_1000>(startPoint, endPoint, spanRadians, 0.0));

        Arc_T<double> arc;
        // convert both endpoints to double precision for each number in the matrix using a for loop
        for (int i = 0; i < 3; i++) {
            arc(i, 0) = startPoint(i).convert_to<double>();
            arc(i, 1) = endPoint(i).convert_to<double>();
        }

        //generate the constant Zs, which is randomly generated between the two endpoints
        double randZ = Dist(gen);
        mpfr_float_1000 constZ_mpfr = startPoint(2) + randZ * abs(endPoint(2) - startPoint(2));
        double constZ = static_cast<double>(constZ_mpfr);
        constZs.push_back(constZ);

        arcs.push_back(arc);
    }

    return arcs;
}



#endif //CODES_ARCS_GENERATOR_H
