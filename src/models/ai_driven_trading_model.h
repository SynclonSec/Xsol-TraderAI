
#pragma once

#include <vector>
#include <memory>
#include "../utils/fft.h"
#include "../utils/stochastic_calculus.h"
#include "../utils/math_utils.h"

namespace xsol {

class AIDrivenTradingModel {
public:
    struct PredictionResult {
        double predictedPrice;
        double confidence;
        double volatility;
        std::vector<double> probabilityDistribution;
    };

    struct ModelConfig {
        size_t windowSize;
        double learningRate;
        double momentumFactor;
        double volatilityThreshold;
        double confidenceThreshold;
    };

    AIDrivenTradingModel(const ModelConfig& config = ModelConfig());
    ~AIDrivenTradingModel();

    // Train the model with historical data
    bool train(const std::vector<double>& prices,
              const std::vector<double>& volumes,
              const std::vector<double>& additionalFeatures);

    // Make predictions
    PredictionResult predict(const std::vector<double>& currentMarketData);

    // Analyze volatility using FFT
    double analyzeVolatility(const std::vector<double>& prices);

    // Update model with new market data
    void updateModel(const std::vector<double>& newData);

private:
    struct Impl;
    std::unique_ptr<Impl> pImpl;
    
    math::FFTAnalyzer m_fftAnalyzer;
    ModelConfig m_config;
};

} // namespace xsol