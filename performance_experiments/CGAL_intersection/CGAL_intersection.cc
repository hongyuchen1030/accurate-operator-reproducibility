#include <Eigen/Dense>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/iterator.h>
#include <tbb/parallel_for.h>

#include "../generate_input.hh"

#include <tbb/parallel_for.h>
#include <tbb/global_control.h>
std::unique_ptr<tbb::global_control> g_global_control;
void set_max_num_tbb_threads(int num_threads) {
    if (num_threads < 1) throw std::runtime_error("num_threads must be >= 1");
    g_global_control = std::make_unique<tbb::global_control>(tbb::global_control::parameter::max_allowed_parallelism, num_threads);
}
 
using SK = CGAL::Exact_spherical_kernel_3;
using PT = typename SK::Point_3;

PT to_cgal(const Eigen::Ref<const Eigen::Vector3d> &p) { return PT(p[0], p[1], p[2]); }

// Convert a CGAL point type into Eigen (rounding each of its components
// to floating point). Note that in the case of an intersection point,
// the input components are algebraic numbers of degree 2.
template<typename CGALPt3>
Eigen::Vector3d to_eigen(const CGALPt3 &p) {
    Eigen::Vector3d result;

    result << CGAL::to_double(p.x()),
              CGAL::to_double(p.y()),
              CGAL::to_double(p.z());

    return result;
}

template<bool Parallel = false>
void compute_values_CGAL(const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    const size_t size = ptsA.rows();

    if ((size != ptsB.rows()) || (size != lattitudes.rows()) || (size != X.rows()) || (size != Y.rows())) {
        throw std::runtime_error("All passed arrays must have the same size");
    }

    PT origin(0, 0, 0);
    SK::Vector_3 z_axis(0, 0, 1);
    SK::Sphere_3 unit_sphere(origin, 1);

    // unsigned integer indicates multiplicity of intersection point
    using Point_and_multiplicity = std::pair<SK::Circular_arc_point_3, unsigned>;

    auto processInput = [&](size_t i) {
        // Great circle passing through A and B.
        SK::Circle_3 c1(unit_sphere, SK::Plane_3(origin, to_cgal(ptsA.row(i)), to_cgal(ptsB.row(i))));

        // Line of constant lattiude
        SK::Circle_3 c2(unit_sphere, SK::Plane_3(PT(0, 0, lattitudes[i]), z_axis));

        SK::Intersect_3 inter;
        std::vector<Point_and_multiplicity> intersections;

        // Only recover point-valued intersections (not circles)
        inter(c1, c2, CGAL::dispatch_or_drop_output<Point_and_multiplicity>(std::back_inserter(intersections)));
        
        const auto &pt = intersections[1]; // The single intersection compted by the other code is the second one computed by CGAL.
        X[i] = CGAL::to_double(intersections[1].first.x());
        Y[i] = CGAL::to_double(intersections[1].first.y());
    };

    if constexpr (Parallel) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, size),
                [&](const tbb::blocked_range<size_t> r) {
                    for (size_t i = r.begin(); i != r.end(); ++i)
                        processInput(i);
                });
    } else {
        for (size_t i = 0; i < size; ++i) { processInput(i); }
    }
}

template <bool Parallel = false>
double run_benchmark_CGAL(const size_t numTests, const Eigen::MatrixXd &ptsA, const Eigen::MatrixXd &ptsB, const Eigen::VectorXd &lattitudes, Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < numTests; ++i) {
         compute_values_CGAL<Parallel>(ptsA, ptsB, lattitudes, X, Y);
    }

    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;

    return elapsed.count();
}

int main() {
    static constexpr size_t dataSize = 1000000;

    Eigen::MatrixXd ptsA, ptsB;
    Eigen::VectorXd lattitudes;
    generateInput(dataSize, ptsA, ptsB, lattitudes);

    // Storage for computed result.
    Eigen::VectorXd X(dataSize), Y(dataSize);

    compute_values_CGAL</* Parallel = */ false>(ptsA, ptsB, lattitudes, X, Y);

    std::cout.precision(18);
    std::cout << "CGAL result: " << X[0] << ", " << Y[0] << ", " << lattitudes[0] << std::endl;

    constexpr size_t numTests = 10;
    double time = run_benchmark_CGAL</* Parallel = */ false>(numTests, ptsA, ptsB, lattitudes, X, Y);
    std::cout << "Serial CGAL time: " << time /  (numTests * dataSize) << " seconds" << std::endl;

    // Parallel versions
    for (size_t numThreads = 1; numThreads < 32; numThreads *= 2) {
        set_max_num_tbb_threads(numThreads);
        time = run_benchmark_CGAL</* Parallel = */ true>(numTests, ptsA, ptsB, lattitudes, X, Y);
        std::cout << "Parallel CGAL (" << numThreads << " threads) time: " << time /  (numTests * dataSize) << " seconds" << std::endl;
    }

    return 0;
}
