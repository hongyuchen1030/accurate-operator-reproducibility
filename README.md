
This repository contains procedures to reproduce the experiments presented in the paper  
**SIAM Journal on Scientific Computing: Fast and Accurate Intersections on a Sphere** by Hongyu Chen, Paul A. Ullrich, and Julian Panetta.


# Accuracy Experiments
All required commands are stored in `./accuracy_experiments/run_accuracy_benchmark.sh`. Since this code uses `libquadmath`, we recommend building with recent version of GCC. The following instructions build the code with GCC-13.


> **Note:** All commands should be executed within the `./accuracy_experiments` directory.  
> Before proceeding, navigate to the correct directory using:
```bash
cd ./accuracy_experiments
```


## Cmake build 
To build the project using `CMake` and `ninja`, run the following command. Make sure to update `-DCGAL_DIR` to your CGAL directory
```bash
# Create and enter the build directory
mkdir -p build && cd build

# Configure the project with CMake
cmake .. -GNinja \
      -DCMAKE_C_COMPILER=gcc-13 \
      -DCMAKE_CXX_COMPILER=g++-13 \
      -DCMAKE_BUILD_TYPE=Release \
      -DCGAL_DIR=<Your-CGAL-Directory>

# Compile the project
ninja

```

## Generate Data Set
Use the following command to generate both the primary and secondary datasets, if needed. The generated data will be stored in `./accuracy_experiments/generated_arcs`, with each floating-point result represented as a pair of significand and exponent values.  

In the script below, `LATITUDES` defines the latitude intervals for generating the primary dataset, while `"$OFFSETS"` specifies the offset intervals for the secondary experiments

It's highly recommanded to run `./sanitize_arcs` for the newly-generated arcs before running the experiments.

```bash
# Specify the number of arcs for each region you want to generate
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


# Optionally run the generate_arcs and sanitize_arcs  executables 
./generate_arcs "$NUM_ARCS" "$LATITUDES" "$OFFSETS"
echo "generate_arcs has been executed."

./sanitize_arcs "$MPFR_PRECISIONS" "$LATITUDES" "$OFFSETS"
echo "sanitize_arcs has been executed."

```


### Run the benchmark
Below will generate the benchmark for both primary and secondary datasets that's being stored in the directory `./accuracy_experiments/benchmark_results` and also the intermediate results will also be stored in `./accuracy_experiments/intermediate_results`. each flaoting point number results will be store as a pair of significant and exponent.

```bash
# The selected MPFR precisions,
MPFR_PRECISIONS="16,17,18,19,20,22,24,26,28,30,32"

# Run all method benchmark except for CGAL
./run_benchmarks "$MPFR_PRECISIONS" "$LATITUDES"
echo "run_benchmarks has been executed."

# Run CGAL benchmark
./CGAL_accuracy "$LATITUDES"
echo "Complete running of CGAL_accuracy"

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

### Krumm's Implementation
Krumm's method is implemented in Python and analyzed separately from other methods.
The implementation files are:

- `./accuracy_experiments/gcsxsc.py`
- `./accuracy_experiments/geometries.py`
- `./accuracy_experiments/utils.py`

To run the Krumm’s method, execute:

```bash
python3 ./accuracy_experiments/gcsxsc.py
```

Ensure that the variable **`base_path`** in the `main` function points to the correct directories.

For analysis, use the `./accuracy_experiments/krumm_accuracy.nb` notebook. Update the following inputs:

1. **For primary dataset:**
```mathematica
arcFilename = latRangeStr <> "Arcs_Exponent.csv";
benchmarkKrummFilename = latRangeStr <> "_krumm_double_ArcUp.csv";
```
2. **For secondary dataset:**
```mathematica
arcFilename = latRangeStr <> "ArcUp_Arcs_Exponent.csv";
benchmarkKrummFilename = latRangeStr <> "_krumm_double_ArcUp.csv";
```

Refer to the naming convention in `./generated_arcs` and `./benchmark_results` directories to ensure correctness.


# Performance Experiments

This section provides instruction for running the performance experiments
> **Note:** All commands should be executed within the `./performance_experiments` directory.  
> Before proceeding, navigate to the correct directory using:
```bash
cd ./performance_experiments
```

## For Direct Floating Points/ MPFR and Quadruple Precisions

### Cmake build 
To build the project using `CMake` , run the following command.  Since this code uses `libquadmath`, we recommend building with recent version of GCC. The following instructions build the code with GCC-13.
```bash
# Create and enter the build directory
mkdir -p build && cd build

# Configure the project with CMake
cmake .. -GNinja \
      -DCMAKE_C_COMPILER=gcc-13 \
      -DCMAKE_CXX_COMPILER=g++-13 \
      -DCMAKE_BUILD_TYPE=Release \

# Compile the project
make || exit 1
```

### Run Benchmarks
The scalar performance benchmark is stored in `./GCA_ConstLat_Performance_Benchmark_scalar_MPFR_Points.cc`. Before running, adjust the following values inside its `main` function as needed.
The output is stored in `./performance_experiments/output/average_time_scalar_MPFR_points.csv`.
```C++
      static constexpr size_t dataSize = 10000000; // Adjust the datasize values
```

The vectorized and parallelized performance benchmark is stored in `./performance_experiments/GCA_ConstLat_Performance_Benchmark_parallelized.cc`. Similarly, adjust the following values:
```C++
      static constexpr size_t numTest_hardcoded = 100; // Number of test repeat to make sure cache is warmed up
      static constexpr size_t dataSize = 100000000;
```
The output is stored in `./performance_experiments/output/average_time_parallel.csv`.

To run both benchmarks, simply use the following commands in the correct directory:

```bash
# Run the scalar benchmark executable
echo "Running the scalar benchmark..."
./GCA_ConstLat_Performance_Benchmark_scalar_MPFR_Points || exit 1

# Run the parallel benchmark executable
echo "Running the parallel benchmark..."
./GCA_ConstLat_Performance_Benchmark_parallelized || exit 1

```


## For CGAL performance benchmark
Navigate into the  `./performance_experiments/CGAL_intersection` directory

### Cmake build 
To build the project using `CMake` and `ninja`, run the following command. For consistency with other experiments, these instructions use `GCC-13`.
Make sure to update `-DCGAL_DIR` with the path to your CGAL directory.
```bash
# Create and enter the build directory
mkdir -p build && cd build

# Configure the project with CMake
cmake .. -GNinja \
      -DCMAKE_C_COMPILER=gcc-13 \
      -DCMAKE_CXX_COMPILER=g++-13 \
      -DCMAKE_BUILD_TYPE=Release \
      -DCGAL_DIR=<Your-CGAL-Directory>

# Compile the project
ninja

```

### Run Benchmarks
All benchmarks, both scalar and parallelized, are located in `./performance_experiments/CGAL_intersection/CGAL_intersection.cc`.

Similarly, adjust the following value in the main function to your preferred setting:
```C++
      static constexpr size_t dataSize = 10000000;
```

To run the performance use 
```bash
./CGAL_intersection
```


The results will be displayed in the terminal.
