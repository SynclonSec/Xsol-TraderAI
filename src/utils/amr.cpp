#include "amr.h"
#include <algorithm>
#include <stdexcept>
#include <numeric>
#include <iostream>
#include <fstream>

namespace xsol {
namespace math {

AdaptiveMeshRefinement::AdaptiveMeshRefinement(const AMRConfig& config)
    : m_config(config) {
    if (config.maxLevels < 1 || config.errorThreshold <= 0) {
        throw std::invalid_argument("Invalid AMR configuration");
    }
    m_meshLevels.resize(config.maxLevels);
}

std::vector<AdaptiveMeshRefinement::MeshPoint> 
AdaptiveMeshRefinement::refineMesh(
    const std::vector<double>& initialData,
    const std::function<double(double)>& errorEstimator) {

    // Initialize base level
    m_meshLevels[0].resize(initialData.size());
    for (size_t i = 0; i < initialData.size(); ++i) {
        m_meshLevels[0][i] = MeshPoint{
            .value = initialData[i],
            .error = 0.0,
            .needsRefinement = false
        };
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

        // Refine marked regions
        for (size_t i = 0; i < m_meshLevels[level].size(); ++i) {
            if (m_meshLevels[level][i].needsRefinement) {
                size_t endIdx = i + 1;
                while (endIdx < m_meshLevels[level].size() && 
                       m_meshLevels[level][endIdx].needsRefinement) {
                    ++endIdx;
                }
                refineRegion(level, i, endIdx);
                i = endIdx - 1;
            }
        }
    }

    // Combine all levels into final mesh
    std::vector<MeshPoint> refinedMesh;
    for (const auto& level : m_meshLevels) {
        refinedMesh.insert(refinedMesh.end(), level.begin(), level.end());
    }

    return refinedMesh;
}

void AdaptiveMeshRefinement::adaptResources(
    const std::vector<double>& activityMetrics) {

    double totalActivity = std::accumulate(
        activityMetrics.begin(), 
        activityMetrics.end(), 
        0.0
    );

    if (totalActivity <= 0) return;

    // Normalize activity metrics
    std::vector<double> normalizedActivity(activityMetrics.size());
    std::transform(
        activityMetrics.begin(), 
        activityMetrics.end(), 
        normalizedActivity.begin(),
        [totalActivity](double x) { return x / totalActivity; }
    );

    // Adjust refinement levels based on activity
    for (size_t i = 0; i < normalizedActivity.size(); ++i) {
        size_t targetLevel = static_cast<size_t>(
            normalizedActivity[i] * m_config.maxLevels
        );
        if (targetLevel >= m_meshLevels.size()) {
            targetLevel = m_meshLevels.size() - 1;
        }
        // Ensure minimum points per level
        while (m_meshLevels[targetLevel].size() < m_config.minPointsPerLevel) {
            refineRegion(targetLevel - 1, 0, m_meshLevels[targetLevel - 1].size());
        }
    }
}

double AdaptiveMeshRefinement::computeLocalError(
    const MeshPoint& point,
    const std::function<double(double)>& errorEstimator) {
    
    return std::abs(errorEstimator(point.value));
}

void AdaptiveMeshRefinement::refineRegion(
    size_t level, 
    size_t startIdx, 
    size_t endIdx) {
    
    if (level >= m_config.maxLevels - 1) return;

    size_t numNewPoints = static_cast<size_t>(
        (endIdx - startIdx) * m_config.refinementRatio
    );

    std::vector<MeshPoint> refinedPoints(numNewPoints);
    double dx = 1.0 / (numNewPoints - 1);

    for (size_t i = 0; i < numNewPoints; ++i) {
        double t = i * dx;
        size_t baseIdx = startIdx + (i * (endIdx - startIdx)) / numNewPoints;
        
        refinedPoints[i] = MeshPoint{
            .value = m_meshLevels[level][baseIdx].value,
            .error = 0.0,
            .needsRefinement = false
        };
    }

    m_meshLevels[level + 1].insert(
        m_meshLevels[level + 1].end(),
        refinedPoints.begin(),
        refinedPoints.end()
    );
}

void AdaptiveMeshRefinement::logMeshState(const std::string& filename) {
    std::ofstream logFile(filename);
    if (!logFile.is_open()) {
        throw std::runtime_error("Unable to open log file");
    }

    for (size_t level = 0; level < m_meshLevels.size(); ++level) {
        logFile << "Level " << level << ":\n";
        for (const auto& point : m_meshLevels[level]) {
            logFile << "Value: " << point.value << ", Error: " << point.error 
                    << ", Needs Refinement: " << (point.needsRefinement ? "Yes" : "No") << "\n";
        }
        logFile << "\n";
    }

    logFile.close();
}

void AdaptiveMeshRefinement::visualizeMesh() const {
    // Placeholder for mesh visualization logic
    // This could involve generating a plot or other visual representation
    for (size_t level = 0; level < m_meshLevels.size(); ++level) {
        std::cout << "Level " << level << ":\n";
        for (const auto& point : m_meshLevels[level]) {
            std::cout << "Value: " << point.value << ", Error: " << point.error 
                      << ", Needs Refinement: " << (point.needsRefinement ? "Yes" : "No") << "\n";
        }
        std::cout << "\n";
    }
}

} // namespace math
} // namespace xsol
