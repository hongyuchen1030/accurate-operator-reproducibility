# Accuracy Experiments
This repository contains procedures to **reproduce the accuracy experiments** presented in the paper  
**SIAM Journal on Scientific Computing: Accurate Intersection Point Calculation**.

All required commands are stored in `./accuracy_experiments/run_accuracy_benchmark.sh`. Since the project relies on `mppp` support, it is recommended to use `GCC 13` for building.

> **Note:** All commands should be executed within the `./accuracy_experiments` directory.  
> Before proceeding, navigate to the correct directory using:
```bash
cd ./accuracy_experiments
```

## Cmake build
To build the project using CMake, run the following command
```bash
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

## Generate Data Set
Use the following command to generate both the primary and secondary datasets, if needed. The generated data will be stored in `./accuracy_experiments/generated_arcs`, with each floating-point result represented as a pair of significant and exponent values.  

In the script below, `LATITUDES` defines the latitude intervals for generating the primary dataset, while `"$OFFSETS"` specifies the offset intervals for the secondary experiments.s

```bash
# Specify the number of arcs you want to generate
NUM_ARCS="100000"


# Generate the preferred lattitude interval
LATITUDES="0,0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,20,30,40,50,60,70,80,89,89.9"


# Generate 5 points between each specified range around pole area
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


./generate_arcs "$NUM_ARCS" "$LATITUDES" "$OFFSETS"
echo "generate_arcs has been executed."

```


### Run the benchmark
Below will generate the benchmark for both primary and secondary datasets that's being stored in the directory `./accuracy_experiments/benchmark_results` and also the intermediate results will also be stored in `./accuracy_experiments/intermediate_results`. each flaoting point number results will be store as a pair of significant and exponent.
```
# The selected MPFR precisions,
MPFR_PRECISIONS="16,17,18,19,20,22,24,26,28,30,32"

./run_benchmarks "$MPFR_PRECISIONS" "$LATITUDES"
echo "run_benchmarks has been executed."

./intermediate_results "$MPFR_PRECISIONS" "$LATITUDES"
echo "intermediate_results has been executed."
```

### Run Analysis 
Now utilizing the Wolframe mathematica to read the benchmark result and perform the analysis, the results will be stored in `./accuracy_experiments/mathematica_data`

```bash
# Now run the mathematica script
WOLFRAMSCRIPT_PATH="<Your WOLFRAMSCRIPT Path>"

# Run AccuracyAnalysis.m
"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/AccuracyAnalysis.m" -- "$MPFR_PRECISIONS" "$LATITUDES"
echo "AccuracyAnalysis.m has been executed."

"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/AccuracyAnalysisArcUp.m" -- "$MPFR_PRECISIONS" "$OFFSETS"
echo "AccuracyAnalysisArcUp.m has been executed."

"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/AccuracyAnalysisArcUpCGAL.m" -- "$MPFR_PRECISIONS" "$OFFSETS"
echo "AccuracyAnalysisArcUpCGAL.m has been executed."


# # Define query ranges
QUERY_RANGE="{0,1},{2,89},{89,90}"
QUERY_RANGE_ENTIRE_NORTH="{0,90}"
# Run IntermediateAnalysis.m for the entire globe
"$WOLFRAMSCRIPT_PATH" -file "$SRC_DIR/IntermediateAnalysis.m" -- "$MPFR_PRECISIONS" "$LATITUDES" "$QUERY_RANGE_ENTIRE_NORTH"
echo "IntermediateAnalysis.m for entire globe has been executed."
```

