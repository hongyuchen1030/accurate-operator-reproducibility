cmake_minimum_required(VERSION 3.5)
project(CGAL_intersection)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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

find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)

add_executable(CGAL_intersection CGAL_intersection.cc)
target_link_libraries(CGAL_intersection CGAL::CGAL Eigen3::Eigen TBB::tbb)

