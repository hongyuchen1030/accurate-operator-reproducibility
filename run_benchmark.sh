#!/bin/bash

# Define the build directory and source directory
BUILD_DIR="build_gcc13"
SRC_DIR="/home/hyvchen/AccuracyBenchmarkEFT"

# Clean up the old build directory if it exists
if [ -d "$BUILD_DIR" ]; then
  echo "Removing old build directory: $BUILD_DIR"
  rm -rf "$BUILD_DIR"
fi

# Create the new build directory
echo "Creating build directory: $BUILD_DIR"
mkdir "$BUILD_DIR"
cd "$BUILD_DIR" || exit 1

# Run CMake with GCC 13 and Release mode
echo "Configuring the project with CMake..."
cmake -DCMAKE_C_COMPILER=gcc-13 -DCMAKE_CXX_COMPILER=g++-13 -DCMAKE_BUILD_TYPE=Release "$SRC_DIR" || exit 1

make


NUM_ARCS="10000"
MPFR_PRECISIONS="16 17"

# Generate 5 points between each specified range
LATITUDES="0,0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,20,30,40,50,60,70,80,89,89.9"
LATITUDES+=$(echo ",89.99,$(seq -f '%.6f' 89.99 0.0002 89.999 | tail -n +2 | head -n 4 | tr '\n' ',' | sed 's/,$//')")
LATITUDES+=$(echo ",89.999,$(seq -f '%.7f' 89.999 0.00004 89.9999 | tail -n +2 | head -n 4 | tr '\n' ',' | sed 's/,$//')")
LATITUDES+=$(echo ",89.9999,$(seq -f '%.8f' 89.9999 0.000008 89.99999 | tail -n +2 | head -n 4 | tr '\n' ',' | sed 's/,$//')")
LATITUDES+=",89.99999,90"
# Define query ranges
QUERY_RANGE="{0,1},{2,89},{89,90}"
QUERY_RANGE_ENTIRE_NORTH="{0,90}"


# # Optionally run the generated C++ executables
# ./generate_arcs "$NUM_ARCS" "$LATITUDES"
# echo "generate_arcs has been executed."

# ./run_benchmarks "$MPFR_PRECISIONS" "$LATITUDES"
# echo "run_benchmarks has been executed."

./intermediate_results  "$MPFR_PRECISIONS" "$LATITUDES"
echo "intermediate_results has been executed."