#pragma once

#include <vector>
#include <memory>
#include "ai_driven_trading_model.h"
#include "wallet_cluster_analysis_model.h"
#include "../utils/llm_inference.h"

namespace xsol {

class MultiAgentSystemModel {
public:
    struct Agent {
        std::string id;
        double performance;
        double confidence;
        std::vector<double> weights;
        AIDrivenTradingModel::PredictionResult lastPrediction;
    };

    struct AgentConfig {
        size_t numAgents;
        double learningRate;
        double explorationRate;
        double adaptationRate;
        size_t historyWindow;
    };

    struct SystemState {
        std::vector<Agent> agents;
        double systemConfidence;
        double consensusLevel;
        std::vector<double> aggregatedPrediction;
    };

    MultiAgentSystemModel(const AgentConfig& config);
    ~MultiAgentSystemModel();

    // Initialize the multi-agent system
    void initialize(size_t inputDimension);

    // Update agent weights based on performance
    void updateAgents(const std::vector<double>& actualOutcomes);

    // Get system prediction
    std::vector<double> getPrediction(const std::vector<double>& marketData);

    // Get current system state
    SystemState getSystemState() const;

private:
    struct Impl;
    std::unique_ptr<Impl> pImpl;
    AgentConfig m_config;

    // Helper methods
    double calculateConsensus(const std::vector<Agent>& agents);
    void evolveAgents();
    std::vector<double> aggregatePredictions(
        const std::vector<AIDrivenTradingModel::PredictionResult>& predictions
    );
};

} // namespace xsol