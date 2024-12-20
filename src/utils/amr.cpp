#include <vector>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <future> // For parallel refinement
#include <mutex>  // For thread-safe operations

namespace xsol {
namespace math {

// Configuration for Adaptive Mesh Refinement
struct AMRConfig {
    size_t maxLevels;         // Maximum number of refinement levels
    double errorThreshold;    // Threshold for refinement
    double refinementRatio;   // Ratio of refinement for dividing intervals
    bool parallelRefinement;  // Whether to enable parallel refinement
};

// A point in the mesh
struct MeshPoint {
    double value; // Data value at this point
    double error; // Error at this point
    bool needsRefinement; // Whether this point needs refinement
};

class AdaptiveMeshRefinement {
public:
    AdaptiveMeshRefinement(const AMRConfig& config)
        : m_config(config) {
        if (config.maxLevels < 1 || config.errorThreshold <= 0 || config.refinementRatio <= 1.0) {
            throw std::invalid_argument("Invalid AMR configuration");
        }
        m_meshLevels.resize(config.maxLevels);
    }

    std::vector<MeshPoint> refineMesh(
        const std::vector<double>& initialData,
        const std::function<double(double)>& errorEstimator);

private:
    void refineMarkedRegions(size_t level);
    void refineRegion(size_t level, size_t startIdx, size_t endIdx);
    double computeLocalError(const MeshPoint& point, const std::function<double(double)>& errorEstimator);
    double cubicSplineInterpolation(const std::vector<double>& x, const std::vector<double>& y, double xi);

    AMRConfig m_config;
    std::vector<std::vector<MeshPoint>> m_meshLevels;
    std::mutex m_mutex; // For thread-safe operations during parallel refinement
};

std::vector<MeshPoint> AdaptiveMeshRefinement::refineMesh(
    const std::vector<double>& initialData,
    const std::function<double(double)>& errorEstimator) {

    if (initialData.empty()) {
        throw std::invalid_argument("Initial data cannot be empty");
    }

    // Initialize base level
    m_meshLevels[0].resize(initialData.size());
    for (size_t i = 0; i < initialData.size(); ++i) {
        m_meshLevels[0][i] = {initialData[i], 0.0, false};
    }

    // Refine mesh levels
    for (size_t level = 0; level < m_config.maxLevels - 1; ++level) {
        bool needsMoreRefinement = false;

        // Compute errors and mark points for refinement
        for (auto& point : m_meshLevels[level]) {
            point.error = computeLocalError(point, errorEstimator);
            point.needsRefinement = point.error > m_config.errorThreshold;
            needsMoreRefinement |= point.needsRefinement;
        }

        if (!needsMoreRefinement) break;

        // Refine marked regions, potentially in parallel
        if (m_config.parallelRefinement) {
            std::vector<std::future<void>> tasks;
            for (size_t i = 0; i < m_meshLevels[level].size(); ++i) {
                if (m_meshLevels[level][i].needsRefinement) {
                    size_t startIdx = i;
                    while (i < m_meshLevels[level].size() && m_meshLevels[level][i].needsRefinement) {
                        ++i;
                    }
                    tasks.emplace_back(std::async(std::launch::async, [this, level, startIdx, i]() {
                        refineRegion(level, startIdx, i);
                    }));
                }
            }
            for (auto& task : tasks) {
                task.get();
            }
        } else {
            refineMarkedRegions(level);
        }
    }

    // Combine all levels into the final mesh
    std::vector<MeshPoint> refinedMesh;
    for (const auto& level : m_meshLevels) {
        refinedMesh.insert(refinedMesh.end(), level.begin(), level.end());
    }
    return refinedMesh;
}

void AdaptiveMeshRefinement::refineMarkedRegions(size_t level) {
    if (level >= m_config.maxLevels - 1) return;

    auto& currentLevel = m_meshLevels[level];
    auto& nextLevel = m_meshLevels[level + 1];

    for (size_t i = 0; i < currentLevel.size(); ++i) {
        if (currentLevel[i].needsRefinement) {
            size_t startIdx = i;
            while (i < currentLevel.size() && currentLevel[i].needsRefinement) {
                ++i;
            }
            refineRegion(level, startIdx, i);
        }
    }
}

void AdaptiveMeshRefinement::refineRegion(size_t level, size_t startIdx, size_t endIdx) {
    if (level >= m_config.maxLevels - 1) return;

    auto& currentLevel = m_meshLevels[level];
    auto& nextLevel = m_meshLevels[level + 1];

    size_t numNewPoints = std::max(
        static_cast<size_t>((endIdx - startIdx) * m_config.refinementRatio), size_t(2));
    double dx = 1.0 / (numNewPoints - 1);

    std::vector<double> x(endIdx - startIdx);
    std::vector<double> y(endIdx - startIdx);
    for (size_t i = startIdx; i < endIdx; ++i) {
        x[i - startIdx] = static_cast<double>(i);
        y[i - startIdx] = currentLevel[i].value;
    }

    // Ensure the next level is properly resized
    nextLevel.reserve(nextLevel.size() + numNewPoints);

    for (size_t i = 0; i < numNewPoints; ++i) {
        double t = i * dx;
        double xi = startIdx + t * (endIdx - startIdx);
        double interpolatedValue = cubicSplineInterpolation(x, y, xi);
        std::lock_guard<std::mutex> lock(m_mutex);
        nextLevel.push_back({interpolatedValue, 0.0, false});
    }
}

double AdaptiveMeshRefinement::computeLocalError(
    const MeshPoint& point,
    const std::function<double(double)>& errorEstimator) {
    return std::abs(errorEstimator(point.value));
}

double AdaptiveMeshRefinement::cubicSplineInterpolation(const std::vector<double>& x, const std::vector<double>& y, double xi) {
    size_t n = x.size();
    if (n < 2) {
        throw std::invalid_argument("At least two points are required for cubic spline interpolation");
    }

    // Find the interval [x[i], x[i+1]] that contains xi
    size_t i = 0;
    while (i < n - 1 && x[i] < xi) {
        ++i;
    }

    // Ensure xi is within the bounds of x
    if (i == 0) i = 1;
    if (i == n) i = n - 1;

    double x0 = x[i - 1];
    double x1 = x[i];
    double y0 = y[i - 1];
    double y1 = y[i];

    double h = x1 - x0;
    double a = y0;
    double b = (y1 - y0) / h;
    double c = (3.0 * (y1 - y0) / (h * h)) - ((2.0 * b) / h);
    double d = ((y1 - y0) / (h * h * h)) - (b / (h * h));

    double dx = xi - x0;
    return a + b * dx + c * dx * dx + d * dx * dx * dx;
}

} // namespace math
} // namespace xsol

// Example usage
int main() {
    xsol::math::AMRConfig config = {3, 0.1, 2.0, true}; // Enable parallel refinement
    xsol::math::AdaptiveMeshRefinement amr(config);

    std::vector<double> initialData = {0.0, 1.0, 2.0, 3.0, 4.0};
    auto errorEstimator = [](double value) { return std::fabs(value - 2.0); }; // Example error function

    try {
        auto refinedMesh = amr.refineMesh(initialData, errorEstimator);

        for (const auto& point : refinedMesh) {
            std::cout << "Value: " << point.value << ", Error: " << point.error
                      << ", Needs Refinement: " << (point.needsRefinement ? "Yes" : "No") << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
