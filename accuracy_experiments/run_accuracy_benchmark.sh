# #!/bin/bash

# Define the build directory and source directory
BUILD_DIR="build_gcc13"
SRC_DIR="/home/hyvchen/accurate-operator-reproducibility/accuracy_experiments"

Clean up the old build directory if it exists
if [ -d "$BUILD_DIR" ]; then
  echo "Removing old build directory: $BUILD_DIR"
  rm -rf "$BUILD_DIR"
fi

# Create the new build directory
echo "Creating build directory: $BUILD_DIR"
mkdir "$BUILD_DIR"
cd "$BUILD_DIR" || exit 1

echo "Configuring the project with CMake..."
cmake -GNinja \
      -DCMAKE_C_COMPILER=gcc-13 \
      -DCMAKE_CXX_COMPILER=g++-13 \
      -DCMAKE_BUILD_TYPE=Release \
      -DCGAL_DIR=/home/jpanetta/software/CGAL-6.0 \
      "$SRC_DIR"
ninja



NUM_ARCS="100"
MPFR_PRECISIONS="16,17,18,19,20,22,24,26,28,30,32"

# Generate 5 points between each specified range
LATITUDES="0,0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,20,30,40,50,60,70,80,89,89.9"
LATITUDES+=$(echo ",89.99,$(seq -f '%.6f' 89.99 0.0002 89.999 | tail -n +2 | head -n 4 | tr '\n' ',' | sed 's/,$//')")
LATITUDES+=$(echo ",89.999,$(seq -f '%.7f' 89.999 0.00004 89.9999 | tail -n +2 | head -n 4 | tr '\n' ',' | sed 's/,$//')")
LATITUDES+=$(echo ",89.9999,$(seq -f '%.8f' 89.9999 0.000008 89.99999 | tail -n +2 | head -n 4 | tr '\n' ',' | sed 's/,$//')")
LATITUDES+=",89.99999,90"

# Initialize the OFFSETS
OFFSETS="0.01"

# Generate decreasing values from 0.01 to 0.0000001 (1e-15), dividing by 10 at each step
value=0.01
while (( $(awk "BEGIN {print ($value >= 0.00000000000001)}") )); do
    value=$(awk "BEGIN {print ($value / 10)}")
    OFFSETS+=",$(printf "%.16f" $value)"  # Format value to 8 decimal places
done


# Print the OFFSETS
echo "Generated OFFSETS: $OFFSETS"


# # Define query ranges
QUERY_RANGE="{0,1},{2,89},{89,90}"
QUERY_RANGE_ENTIRE_NORTH="{0,90}"


# Optionally run the generate_arcs and sanitize_arcs  executables 
# ./generate_arcs "$NUM_ARCS" "$LATITUDES" "$OFFSETS"
# echo "generate_arcs has been executed."

# ./sanitize_arcs "$MPFR_PRECISIONS" "$LATITUDES" "$OFFSETS"
# echo "sanitize_arcs has been executed."

./run_benchmarks "$MPFR_PRECISIONS" "$LATITUDES" "$OFFSETS"
echo "run_benchmarks has been executed."

./intermediate_results "$MPFR_PRECISIONS" "$LATITUDES"
echo "intermediate_results has been executed."


# Run CGAL benchmark
./CGAL_accuracy "$LATITUDES" "$OFFSETS"
echo "Complete running of CGAL_accuracy"

# Now run the mathematica script
WOLFRAMSCRIPT_PATH="/usr/bin/wolframscript"

# Run AccuracyAnalysis.m for all methods except for CGAL
"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/AccuracyAnalysis.m" -- "$MPFR_PRECISIONS" "$LATITUDES"
echo "AccuracyAnalysis.m has been executed."

# Run CGAL analysis
"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/AccuracyAnalysisCGAL.m"  -- "$LATITUDES"
echo "AccuracyAnalysisCGAL.m has been executed."

# Run IntermediateAnalysis.m for different locations
"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/IntermediateAnalysis.m" -- "$MPFR_PRECISIONS" "$LATITUDES" "$QUERY_RANGE"
echo "IntermediateAnalysis.m for different locations has been executed."

# Run IntermediateAnalysis.m for the entire globe
"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/IntermediateAnalysis.m" -- "$MPFR_PRECISIONS" "$LATITUDES" "$QUERY_RANGE_ENTIRE_NORTH"
echo "IntermediateAnalysis.m for entire globe has been executed."

"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/AccuracyAnalysisArcUp.m" -- "$MPFR_PRECISIONS" "$OFFSETS"
echo "AccuracyAnalysisArcUp.m has been executed."

"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/AccuracyAnalysisArcUpCGAL.m" -- "$MPFR_PRECISIONS" "$OFFSETS"
echo "AccuracyAnalysisArcUpCGAL.m has been executed."

