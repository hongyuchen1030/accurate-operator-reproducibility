# VectorizedEFT
Utilize Eigen Vectorization for EFT Methods


## Project Structure

- **GCA_ConstLat_Performance_Benchmark.cpp**: Contains the main program logic, where different `VEC_WIDTH` size Eigen arrays are constructed as inputs for the benchmark. Notably, when `vec_width = 1`, it runs the scalar version:

    ```cpp
    for (int i = 0; i < iteration; i++) {
        auto res = gcaHelper.gca_constLat_intersection_accurate_performance(pointA_arr, pointB_arr, constZ);
        USE_RESULT(res);
    }
    ```

- **GCAConstLatIntersections.h**: Contains the core intersection point calculation logic. This header implements the intersection and EFT methods using template-based programming to handle both scalar and vectorized operations.

- **simd_fma.hh**: A file that sets up SIMD Fused Multiply-Add (FMA) operations for the performance testing.

## How to Run

### 1. Running the Benchmark

You can directly run the benchmark using the script `./run_analysis.sh`. 

Inside `run_analysis.sh`, you can specify the following settings to control the benchmark:

```
# Define the iteration count
ITERATION="10"

# Define the vectorization width (vec_width)
VEC_WIDTH="2"

# Define MPFR precisions
MPFR_PRECISIONS="16 17"
```


### 2. Reproduce the Problem
To reproduce the issue of vectorization performance, simply modify the vec_width to be greater than 1 in the `run_analysis.sh` script. The terminal output will display the average execution time for both the floating-point equation implementation and the EFT method.

By setting `VEC_WIDTH` > 1, you can compare the performance of scalar versus vectorized operations.

### 3. Output

The program prints out the average time for both the floating-point implementation and the EFT method in the terminal. You can also navigate to the `./benchmark_results/` directory to review the saved results in `performance_results.csv`.

