# AccuracyBenchmarkEFT
This repository contains the procedure to reproduce the accuracy experiment in the paper SIAM Journal on Scientific Computing: Accurate Intersection Point Calculation  

All the required command are stored in the `./run_benchmark.sh`, due to the `mppp` support, it's recommand to use `gcc13` when building the project.

## Cmake build
To build the project using CMake, run the following command
```
BUILD_DIR="build_gcc13"
SRC_DIR="<source_directory>"


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

```

## Primary Data Set

### Generate Primary Data Set
Use the following command to generate the primary dataset if needed and being stored in the directory `./generated_arcs`. each flaoting point number results will be store as a pair of significant and exponent
```
# Specify the number of arcs you want to generate
NUM_ARCS="100000"


# Generate the preferred lattitude interval
LATITUDES="0,0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,20,30,40,50,60,70,80,89,89.9"


# Generate 5 points between each specified range around pole area
LATITUDES+=$(echo ",89.99,$(seq -f '%.6f' 89.99 0.0002 89.999 | tail -n +2 | head -n 4 | tr '\n' ',' | sed 's/,$//')")
LATITUDES+=$(echo ",89.999,$(seq -f '%.7f' 89.999 0.00004 89.9999 | tail -n +2 | head -n 4 | tr '\n' ',' | sed 's/,$//')")
LATITUDES+=$(echo ",89.9999,$(seq -f '%.8f' 89.9999 0.000008 89.99999 | tail -n +2 | head -n 4 | tr '\n' ',' | sed 's/,$//')")
LATITUDES+=",89.99999,90"

# # Define query ranges
QUERY_RANGE="{0,1},{2,89},{89,90}"
QUERY_RANGE_ENTIRE_NORTH="{0,90}"

# ./generate_arcs "$NUM_ARCS" "$LATITUDES"
# echo "generate_arcs has been executed."

```

### Run the primary dataset benchmark
Below will generate the benchmark results that's being stored in the directory `./benchmark_results` and also the intermediate results will also be stored in `./intermediate_results`. each flaoting point number results will be store as a pair of significant and exponent.
```
# The selected MPFR precisions,
MPFR_PRECISIONS="16,17,18,19,20,22,24,26,28,30,32"

./run_benchmarks "$MPFR_PRECISIONS" "$LATITUDES"
echo "run_benchmarks has been executed."

./intermediate_results "$MPFR_PRECISIONS" "$LATITUDES"
echo "intermediate_results has been executed."
```

### Run Analysis on the primary dataset
Now utilizing the Wolframe mathematica to read the benchmark result and perform the analysis, the results will be stored in `mathematica_data`

```
# Now run the mathematica script
WOLFRAMSCRIPT_PATH="<Your WOLFRAMSCRIPT Path>"

# Run AccuracyAnalysis.m
"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/AccuracyAnalysis.m" -- "$MPFR_PRECISIONS" "$LATITUDES"
echo "AccuracyAnalysis.m has been executed."

# # Define query ranges
QUERY_RANGE="{0,1},{2,89},{89,90}"
QUERY_RANGE_ENTIRE_NORTH="{0,90}"
# Run IntermediateAnalysis.m for the entire globe
"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/IntermediateAnalysis.m" -- "$MPFR_PRECISIONS" "$LATITUDES" "$QUERY_RANGE_ENTIRE_NORTH"
echo "IntermediateAnalysis.m for entire globe has been executed."
```

## Generate Secondary Data Set (the ill-condition case)
Use the following command to generate the secondary sdataset if needed
```
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
```