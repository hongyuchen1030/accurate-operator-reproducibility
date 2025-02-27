cmake_minimum_required(VERSION 3.5)
project(AccuracyBenchmarkEFT)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Add custom CMake module path for FindMPFR.cmake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Find Eigen3
find_package(Eigen3 REQUIRED)

# Find MPFR, GMP, and MP++ libraries
find_package(MPFR REQUIRED)
find_package(GMP REQUIRED)
find_package(mp++ REQUIRED)
find_package(CGAL REQUIRED)

# Create simd_fma library
add_library(simd_fma INTERFACE)
target_link_libraries(simd_fma INTERFACE Eigen3::Eigen)

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


add_executable(generate_arcs generate_arcs.cpp)
target_link_libraries(generate_arcs
    simd_fma
    MPFR::MPFR
    GMP::GMP
    mp++::mp++)


add_executable(run_benchmarks run_benchmarks.cpp)
target_link_libraries(run_benchmarks
    simd_fma
    MPFR::MPFR
    GMP::GMP
    mp++::mp++)

add_executable(intermediate_results intermediate_results.cpp)
target_link_libraries(intermediate_results
    simd_fma
    MPFR::MPFR
    GMP::GMP
    mp++::mp++)


add_executable(CGAL_accuracy CGAL_accuracy.cc)
target_link_libraries(CGAL_accuracy CGAL::CGAL Eigen3::Eigen)
        
    