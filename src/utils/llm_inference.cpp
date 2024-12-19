
#include "llm_inference.h"
#include <curl/curl.h>
#include <nlohmann/json.hpp>
#include <sstream>
#include <regex>

namespace xsol {
namespace ai {

struct LLMInference::Impl {
    CURL* curl;
    std::string lastError;
    std::vector<InferenceResult> cache;
    
    static size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
        ((std::string*)userp)->append((char*)contents, size * nmemb);
        return size * nmemb;
    }
};

LLMInference::LLMInference(const ModelConfig& config) 
    : m_config(config), pImpl(std::make_unique<Impl>()) {
    
    pImpl->curl = curl_easy_init();
    if (!pImpl->curl) {
        throw std::runtime_error("Failed to initialize CURL");
    }
}

LLMInference::~LLMInference() {
    if (pImpl->curl) {
        curl_easy_cleanup(pImpl->curl);
    }
}

std::string LLMInference::formatPrompt(const MarketContext& context) {
    std::stringstream ss;
    ss << "Analyze the following market context:\n\n";
    ss << "Market: " << context.market << "\n";
    ss << "Current Price: " << context.currentPrice << "\n";
    ss << "24h Volume: " << context.volume24h << "\n";
    
    ss << "\nRecent Price History:\n";
    for (const auto& price : context.recentPrices) {
        ss << price << ", ";
    }
    
    ss << "\n\nRecent News:\n";
    for (const auto& news : context.recentNews) {
        ss << "- " << news << "\n";
    }
    
    ss << "\nMarket Sentiment: " << context.marketSentiment << "\n";
    ss << "\nProvide trading analysis and predictions based on this data.";
    
    return ss.str();
}

double LLMInference::calculateConfidence(const std::string& response) {
    // Basic confidence scoring based on response characteristics
    double confidence = 0.5; // Base confidence
    
    // Check for quantitative analysis
    if (response.find("probability") != std::string::npos ||
        response.find("percentage") != std::string::npos) {
        confidence += 0.1;
    }
    
    // Check for reasoning
    if (response.find("because") != std::string::npos ||
        response.find("due to") != std::string::npos) {
        confidence += 0.1;
    }
    
    // Check for market indicators
    if (response.find("indicator") != std::string::npos ||
        response.find("pattern") != std::string::npos) {
        confidence += 0.1;
    }
    
    // Check for multiple factors
    if (response.find("furthermore") != std::string::npos ||
        response.find("additionally") != std::string::npos) {
        confidence += 0.1;
    }
    
    // Cap confidence at 0.95
    return std::min(0.95, confidence);
}

LLMInference::InferenceResult LLMInference::analyzeTradingContext(
    const MarketContext& context) {
    
    InferenceResult result;
    std::string prompt = formatPrompt(context);
    std::string response;
    
    // Prepare API request
    struct curl_slist* headers = nullptr;
    headers = curl_slist_append(headers, "Content-Type: application/json");
    headers = curl_slist_append(headers, 
        ("Authorization: Bearer " + m_config.apiKey).c_str());
    
    nlohmann::json requestBody = {
        {"model", m_config.allowedModels[0]},
        {"prompt", prompt},
        {"max_tokens", 1000},
        {"temperature", 0.7},
        {"top_p", 0.95},
        {"frequency_penalty", 0.0},
        {"presence_penalty", 0.0}
    };

    std::string requestStr = requestBody.dump();
    std::string responseStr;

    curl_easy_setopt(pImpl->curl, CURLOPT_URL, m_config.apiEndpoint.c_str());
    curl_easy_setopt(pImpl->curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(pImpl->curl, CURLOPT_POSTFIELDS, requestStr.c_str());
    curl_easy_setopt(pImpl->curl, CURLOPT_WRITEFUNCTION, Impl::WriteCallback);
    curl_easy_setopt(pImpl->curl, CURLOPT_WRITEDATA, &responseStr);

    CURLcode res = curl_easy_perform(pImpl->curl);
    curl_slist_free_all(headers);

    if (res != CURLE_OK) {
        pImpl->lastError = curl_easy_strerror(res);
        return InferenceResult{false, "", 0.0, "API request failed: " + pImpl->lastError};
    }

    try {
        nlohmann::json responseJson = nlohmann::json::parse(responseStr);
        
        // Extract the model's response
        std::string modelResponse = responseJson["choices"][0]["text"];
        
        // Process the response
        result.success = true;
        result.response = modelResponse;
        result.confidence = calculateConfidence(modelResponse);
        
        // Extract trading signals
        result.signals = extractTradingSignals(modelResponse);
        
        // Cache the result
        if (pImpl->cache.size() >= m_config.maxCacheSize) {
            pImpl->cache.erase(pImpl->cache.begin());
        }
        pImpl->cache.push_back(result);
        
    } catch (const std::exception& e) {
        result.success = false;
        result.error = "Failed to parse API response: " + std::string(e.what());
    }

    return result;
}

std::vector<LLMInference::TradingSignal> 
LLMInference::extractTradingSignals(const std::string& response) {
    std::vector<TradingSignal> signals;
    
    // Regular expressions for pattern matching
    std::regex buyPattern("(?i)\\b(buy|long|bullish)\\b.*?\\$?(\\d+(\\.\\d+)?)");
    std::regex sellPattern("(?i)\\b(sell|short|bearish)\\b.*?\\$?(\\d+(\\.\\d+)?)");
    std::regex confidencePattern("(?i)(confidence|probability).*?(\\d+(\\.\\d+)?)%?");
    
    // Extract buy signals
    for (std::sregex_iterator i(response.begin(), response.end(), buyPattern);
         i != std::sregex_iterator(); ++i) {
        TradingSignal signal;
        signal.action = TradingSignal::Action::BUY;
        signal.price = std::stod((*i)[2]);
        signal.confidence = 0.7; // Default confidence
        signals.push_back(signal);
    }
    
    // Extract sell signals
    for (std::sregex_iterator i(response.begin(), response.end(), sellPattern);
         i != std::sregex_iterator(); ++i) {
        TradingSignal signal;
        signal.action = TradingSignal::Action::SELL;
        signal.price = std::stod((*i)[2]);
        signal.confidence = 0.7; // Default confidence
        signals.push_back(signal);
    }
    
    // Update confidence levels if explicitly mentioned
    std::smatch confidenceMatch;
    if (std::regex_search(response, confidenceMatch, confidencePattern)) {
        double confidence = std::stod(confidenceMatch[2]) / 100.0;
        for (auto& signal : signals) {
            signal.confidence = confidence;
        }
    }
    
    return signals;
}

bool LLMInference::validateResponse(const std::string& response) {
    if (response.empty()) {
        return false;
    }
    
    // Check for minimum response length
    if (response.length() < 50) {
        return false;
    }
    
    // Check for presence of key analysis components
    bool hasAnalysis = response.find("analysis") != std::string::npos ||
                      response.find("indicator") != std::string::npos;
    bool hasPrediction = response.find("predict") != std::string::npos ||
                        response.find("expect") != std::string::npos;
    bool hasReasoning = response.find("because") != std::string::npos ||
                       response.find("due to") != std::string::npos;
    
    return hasAnalysis && hasPrediction && hasReasoning;
}

void LLMInference::updateConfig(const ModelConfig& newConfig) {
    m_config = newConfig;
}

std::string LLMInference::getLastError() const {
    return pImpl->lastError;
}

const std::vector<LLMInference::InferenceResult>& 
LLMInference::getCache() const {
    return pImpl->cache;
}

} // namespace ai
} // namespace xsol