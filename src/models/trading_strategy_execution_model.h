
#pragma once

#include <vector>
#include <string>
#include <memory>
#include "../solana/solana_framework.h"

namespace xsol {

class TradingStrategyExecutionModel {
public:
    struct TradeSignal {
        enum class Action { BUY, SELL, HOLD };
        Action action;
        double quantity;
        double price;
        double confidence;
        std::string market;
        std::string rationale;
    };

    struct ExecutionResult {
        bool success;
        std::string transactionId;
        double executedPrice;
        double executedQuantity;
        double fees;
        std::string status;
        std::string error;
    };

    struct StrategyConfig {
        double maxPositionSize;
        double maxDrawdown;
        double stopLossPercent;
        double takeProfitPercent;
        int maxOpenPositions;
        double minConfidence;
    };

    TradingStrategyExecutionModel(const StrategyConfig& config);
    ~TradingStrategyExecutionModel();

    ExecutionResult executeStrategy(const TradeSignal& signal);
    bool validateSignal(const TradeSignal& signal);
    void updateRiskParameters(const StrategyConfig& newConfig);
    std::vector<ExecutionResult> getOpenPositions();

private:
    struct Impl;
    std::unique_ptr<Impl> pImpl;
    std::shared_ptr<SolanaFramework> m_solanaFramework;
    StrategyConfig m_config;

    bool checkRiskLimits(const TradeSignal& signal);
    double calculateOptimalPositionSize(const TradeSignal& signal);
};

} // namespace xsol