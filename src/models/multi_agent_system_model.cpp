

#include "multi_agent_system_model.h"
#include <random>
#include <algorithm>
#include <cmath>

namespace xsol {

struct MultiAgentSystemModel::Impl {
    std::vector<Agent> agents;
    std::mt19937 rng;
    std::vector<std::vector<double>> predictionHistory;
    size_t inputDimension;
    double systemPerformance;
    
    Impl() : rng(std::random_device{}()) {}
};

MultiAgentSystemModel::MultiAgentSystemModel(const AgentConfig& config)
    : m_config(config), pImpl(std::make_unique<Impl>()) {}

MultiAgentSystemModel::~MultiAgentSystemModel() = default;

void MultiAgentSystemModel::initialize(size_t inputDimension) {
    pImpl->inputDimension = inputDimension;
    pImpl->agents.clear();
    pImpl->agents.reserve(m_config.numAgents);
    
    std::uniform_real_distribution<double> dist(-0.1, 0.1);
    
    // Initialize agents with random weights
    for (size_t i = 0; i < m_config.numAgents; ++i) {
        Agent agent;
        agent.id = "Agent_" + std::to_string(i);
        agent.performance = 0.0;
        agent.confidence = 1.0 / m_config.numAgents;
        agent.weights.resize(inputDimension);
        
        // Initialize weights with small random values
        for (size_t j = 0; j < inputDimension; ++j) {
            agent.weights[j] = dist(pImpl->rng);
        }
        
        pImpl->agents.push_back(agent);
    }
}

double MultiAgentSystemModel::calculateConsensus(
    const std::vector<Agent>& agents) {
    
    if (agents.empty()) return 0.0;
    
    double consensusSum = 0.0;
    int consensusPairs = 0;
    
    // Calculate average agreement between agent predictions
    for (size_t i = 0; i < agents.size(); ++i) {
        for (size_t j = i + 1; j < agents.size(); ++j) {
            double agreement = 0.0;
            const auto& pred1 = agents[i].lastPrediction;
            const auto& pred2 = agents[j].lastPrediction;
            
            // Compare predictions using cosine similarity
            double dotProduct = 0.0;
            double norm1 = 0.0;
            double norm2 = 0.0;
            
            for (size_t k = 0; k < pred1.probabilityDistribution.size(); ++k) {
                dotProduct += pred1.probabilityDistribution[k] * 
                            pred2.probabilityDistribution[k];
                norm1 += pred1.probabilityDistribution[k] * 
                        pred1.probabilityDistribution[k];
                norm2 += pred2.probabilityDistribution[k] * 
                        pred2.probabilityDistribution[k];
            }
            
            if (norm1 > 0 && norm2 > 0) {
                agreement = dotProduct / (std::sqrt(norm1) * std::sqrt(norm2));
            }
            
            consensusSum += agreement;
            consensusPairs++;
        }
    }
    
    return consensusPairs > 0 ? consensusSum / consensusPairs : 0.0;
}

void MultiAgentSystemModel::evolveAgents() {
    std::vector<Agent> newAgents;
    newAgents.reserve(m_config.numAgents);
    
    // Sort agents by performance
    std::sort(pImpl->agents.begin(), pImpl->agents.end(),
        [](const Agent& a, const Agent& b) {
            return a.performance > b.performance;
        }
    );
    
    // Keep top performing agents
    size_t eliteCount = m_config.numAgents / 4;
    for (size_t i = 0; i < eliteCount; ++i) {
        newAgents.push_back(pImpl->agents[i]);
    }
    
    // Create new agents through crossover and mutation
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::normal_distribution<double> mutation(0.0, 0.1);
    
    while (newAgents.size() < m_config.numAgents) {
        // Select parents using tournament selection
        size_t parent1Idx = rand() % eliteCount;
        size_t parent2Idx = rand() % eliteCount;
        
        Agent child = pImpl->agents[parent1Idx];
        
        // Crossover
        for (size_t i = 0; i < child.weights.size(); ++i) {
            if (dist(pImpl->rng) < 0.5) {
                child.weights[i] = pImpl->agents[parent2Idx].weights[i];
            }
            
            // Mutation
            if (dist(pImpl->rng) < m_config.explorationRate) {
                child.weights[i] += mutation(pImpl->rng);
            }
        }
        
        child.performance = 0.0;
        child.confidence = 1.0 / m_config.numAgents;
        newAgents.push_back(child);
    }
    
    pImpl->agents = std::move(newAgents);
}

std::vector<double> MultiAgentSystemModel::aggregatePredictions(
    const std::vector<AIDrivenTradingModel::PredictionResult>& predictions) {
    
    if (predictions.empty()) return {};
    
    size_t predSize = predictions[0].probabilityDistribution.size();
    std::vector<double> aggregated(predSize, 0.0);
    double totalConfidence = 0.0;
    
    // Weighted average based on agent confidence
    for (size_t i = 0; i < pImpl->agents.size(); ++i) {
        const auto& pred = predictions[i];
        double weight = pImpl->agents[i].confidence;
        
        for (size_t j = 0; j < predSize; ++j) {
            aggregated[j] += weight * pred.probabilityDistribution[j];
        }
        totalConfidence += weight;
    }
    
    // Normalize
    if (totalConfidence > 0) {
        for (auto& val : aggregated) {
            val /= totalConfidence;
        }
    }
    
    return aggregated;
}

void MultiAgentSystemModel::updateAgents(
    const std::vector<double>& actualOutcomes) {
    
    // Update agent performance based on prediction accuracy
    for (auto& agent : pImpl->agents) {
        double error = 0.0;
        for (size_t i = 0; i < actualOutcomes.size(); ++i) {
            double pred = agent.lastPrediction.predictedPrice;
            error += std::pow(pred - actualOutcomes[i], 2);
        }
        error = std::sqrt(error / actualOutcomes.size()); // RMSE
        
        // Update performance using exponential moving average
        agent.performance = (1 - m_config.adaptationRate) * agent.performance +
                          m_config.adaptationRate * (1.0 / (1.0 + error));
                          
        // Update confidence
        agent.confidence = agent.performance / 
            std::accumulate(pImpl->agents.begin(), pImpl->agents.end(), 0.0,
                [](double sum, const Agent& a) { return sum + a.performance; });
    }
    
    // Store prediction history
    if (pImpl->predictionHistory.size() >= m_config.historyWindow) {
        pImpl->predictionHistory.erase(pImpl->predictionHistory.begin());
    }
    
    std::vector<double> currentPredictions;
    for (const auto& agent : pImpl->agents) {
        currentPredictions.push_back(agent.lastPrediction.predictedPrice);
    }
    pImpl->predictionHistory.push_back(currentPredictions);
    
    // Evolve agents if needed
    if (pImpl->systemPerformance < 0.5) { // Threshold for evolution
        evolveAgents();
    }
}

std::vector<double> MultiAgentSystemModel::getPrediction(
    const std::vector<double>& marketData) {
    
    std::vector<AIDrivenTradingModel::PredictionResult> predictions;
    predictions.reserve(pImpl->agents.size());
    
    // Get predictions from all agents
    for (auto& agent : pImpl->agents) {
        // Apply agent's weights to market data
        double predictedValue = 0.0;
        for (size_t i = 0; i < marketData.size(); ++i) {
            predictedValue += agent.weights[i] * marketData[i];
        }
        
        AIDrivenTradingModel::PredictionResult result;
        result.predictedPrice = predictedValue;
        result.confidence = agent.confidence;
        
        // Generate probability distribution
        std::normal_distribution<double> dist(predictedValue, 
                                            1.0 / agent.confidence);
        
        result.probabilityDistribution.resize(100);
        for (size_t i = 0; i < 100; ++i) {
            result.probabilityDistribution[i] = dist(pImpl->rng);
        }
        
        agent.lastPrediction = result;
        predictions.push_back(result);
    }
    
    return aggregatePredictions(predictions);
}

MultiAgentSystemModel::SystemState MultiAgentSystemModel::getSystemState() const {
    SystemState state;
    state.agents = pImpl->agents;
    state.consensusLevel = calculateConsensus(pImpl->agents);
    
    // Calculate system confidence as weighted average of agent confidences
    state.systemConfidence = 0.0;
    for (const auto& agent : pImpl->agents) {
        state.systemConfidence += agent.confidence * agent.performance;
    }
    
    // Get latest aggregated prediction
    if (!pImpl->predictionHistory.empty()) {
        state.aggregatedPrediction = pImpl->predictionHistory.back();
    }
    
    return state;
}

} // namespace xsol