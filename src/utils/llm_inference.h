

#pragma once

#include <string>
#include <vector>
#include <memory>
#include <functional>

namespace xsol {
namespace llm {

class LLMInference {
public:
    struct MarketContext {
        std::string market;
        double currentPrice;
        double volume24h;
        std::vector<double> recentPrices;
        std::vector<std::string> recentNews;
        std::vector<std::string> technicalIndicators;
    };

    struct InferenceResult {
        double sentimentScore;
        std::vector<std::string> tradingSignals;
        std::vector<std::pair<std::string, double>> probabilityDistribution;
        std::string reasoning;
        double confidence;
    };

    struct LLMConfig {
        std::string modelEndpoint;
        std::string apiKey;
        size_t maxTokens;
        double temperatureParam;
        double topP;
        int maxRetries;
        double timeoutSeconds;
    };

    LLMInference(const LLMConfig& config);
    ~LLMInference();

    InferenceResult analyzeTradingContext(const MarketContext& context);
    void updateModel(const std::string& newModelVersion);
    double evaluateSentiment(const std::vector<std::string>& texts);

private:
    struct Impl;
    std::unique_ptr<Impl> pImpl;
    LLMConfig m_config;

    std::string formatPrompt(const MarketContext& context);
    InferenceResult parseResponse(const std::string& response);
    bool validateResponse(const InferenceResult& result);
};

} // namespace llm
} // namespace xsol