cmake_minimum_required(VERSION 3.5)
project(VectorizedEFT)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Add custom CMake module path for FindMPFR.cmake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find Eigen3
find_package(Eigen3 REQUIRED)

if (NOT TARGET TBB::tbb)
    include(FetchContent)

    # TBB's tests are causing compilation problems with gcc-13,
    # triggering warnings that are treated as errors. This happens
    # even on the latest versions that purport to fix this.
    #   https://github.com/uxlfoundation/oneTBB/issues/843
    #   https://github.com/uxlfoundation/oneTBB/issues/823
    # So we just disable the tests...
    set(TBB_TEST OFF CACHE BOOL " " FORCE)

    FetchContent_Declare(
        TBB
        GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
        GIT_TAG d75ea937f7e480bd7775d7dab224e0b38451bb40
    )

    FetchContent_MakeAvailable(TBB)
endif()

# Find MPFR, GMP, and MP++ libraries
find_package(MPFR REQUIRED)
find_package(GMP REQUIRED)
find_package(mp++ REQUIRED)

# Create simd_fma library
add_library(simd_fma INTERFACE)
target_link_libraries(simd_fma INTERFACE Eigen3::Eigen)
target_link_libraries(simd_fma INTERFACE TBB::tbb)

# Compilation options for FMA instructions
target_compile_options(simd_fma INTERFACE -march=native -ffp-contract=off)
# A note on -ffp-contract:
#   The default on GCC (at least with -std=gnu17 instead of -std=c++17) is
#   -fp-contract=fast, which can contract floating point multiplications and
#   additions spanning across **multiple statements** into FMA instructions.
#   This can be surprising and change the roundoff behavior. This is especially
#   unsafe for EFT algorithms, which rely fundamentally on the IEEE 754
#   roundoff behavior implied by the sequence of statements as literally
#   written in the source file.



# Add parallized benchmark executable
add_executable(GCA_ConstLat_Performance_Benchmark_parallelized GCA_ConstLat_Performance_Benchmark_parallelized.cc
               instantiations/double_parallel.cc
               instantiations/vec_width2_parallel.cc
               instantiations/vec_width4_parallel.cc
               instantiations/vec_width8_parallel.cc)
target_link_libraries(GCA_ConstLat_Performance_Benchmark_parallelized simd_fma MPFR::MPFR GMP::GMP)

# Add scalar benchmark executable with MPFR Matrix
add_executable(GCA_ConstLat_Performance_Benchmark_scalar_MPFR_Points GCA_ConstLat_Performance_Benchmark_scalar_MPFR_Points.cc)
target_link_libraries(GCA_ConstLat_Performance_Benchmark_scalar_MPFR_Points 
    simd_fma
    MPFR::MPFR
    GMP::GMP
    mp++::mp++)


