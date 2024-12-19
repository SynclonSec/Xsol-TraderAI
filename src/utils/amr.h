#pragma once

#include <vector>
#include <memory>
#include <cmath>
#include <functional>

namespace xsol {
namespace math {

class AdaptiveMeshRefinement {
public:
    struct MeshPoint {
        double value;
        double error;
        bool needsRefinement;
    };

    struct AMRConfig {
        double errorThreshold;
        size_t maxLevels;
        double refinementRatio;
        size_t minPointsPerLevel;
    };

    AdaptiveMeshRefinement(const AMRConfig& config = AMRConfig{
        .errorThreshold = 1e-6,
        .maxLevels = 5,
        .refinementRatio = 2.0,
        .minPointsPerLevel = 10
    });

    // Refine mesh based on error estimators
    std::vector<MeshPoint> refineMesh(
        const std::vector<double>& initialData,
        const std::function<double(double)>& errorEstimator
    );

    // Adapt computational resources based on activity
    void adaptResources(const std::vector<double>& activityMetrics);

private:
    AMRConfig m_config;
    std::vector<std::vector<MeshPoint>> m_meshLevels;

    double computeLocalError(const MeshPoint& point, 
                           const std::function<double(double)>& errorEstimator);
    void refineRegion(size_t level, size_t startIdx, size_t endIdx);
};

} // namespace math
} // namespace xsol