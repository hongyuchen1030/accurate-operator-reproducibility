cmake_minimum_required(VERSION 3.5)
project(CGAL_intersection)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)

FetchContent_Declare(
    TBB
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
    GIT_TAG c9497714821c3d443ee44c732609eb6850195ffb
)

FetchContent_MakeAvailable(TBB)

find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)

add_executable(CGAL_intersection CGAL_intersection.cc)
target_link_libraries(CGAL_intersection CGAL::CGAL Eigen3::Eigen TBB::tbb)

