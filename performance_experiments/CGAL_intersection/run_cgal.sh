#!/bin/bash

# Define the build directory
BUILD_DIR="build_gcc13"
SRC_DIR="/home/hyvchen/accurate-operator-reproducibility/performance_experiments/CGAL_intersection"

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
cmake -GNinja \
      -DCMAKE_C_COMPILER=gcc-13 \
      -DCMAKE_CXX_COMPILER=g++-13 \
      -DCMAKE_BUILD_TYPE=Release \
      -DCGAL_DIR=/home/jpanetta/software/CGAL-6.0 \
      "$SRC_DIR"
ninja


./CGAL_intersection
echo "Complete running of CGAL_intersection"


