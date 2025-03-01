#!/bin/bash

# Define the build directory
BUILD_DIR="build_gcc13"
SRC_DIR="/home/hyvchen/accurate-operator-reproducibility/performance_experiments"

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
cmake -DCMAKE_C_COMPILER=gcc-13 \
      -DCMAKE_CXX_COMPILER=g++-13 \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS="-Wno-error=uninitialized" \
      "$SRC_DIR" || exit 1

# Build the scalar benchmark target
echo "Building the scalar benchmark..."
make || exit 1

# # Run the scalar benchmark executable
# echo "Running the scalar benchmark..."
# ./GCA_ConstLat_Performance_Benchmark_scalar_MPFR_Points || exit 1

# echo "Benchmark run completed!"
