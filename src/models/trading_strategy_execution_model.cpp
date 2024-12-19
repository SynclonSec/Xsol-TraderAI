
#include "trading_strategy_execution_model.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace xsol {

struct TradingStrategyExecutionModel::Impl {
    std::vector<ExecutionResult> openPositions;
    double totalExposure;
    double peakValue;
    double currentDrawdown;
};

TradingStrategyExecutionModel::TradingStrategyExecutionModel(
    const StrategyConfig& config)
    : m_config(config), pImpl(std::make_unique<Impl>()) {
    
    m_solanaFramework = std::make_shared<SolanaFramework>();
    pImpl->totalExposure = 0.0;
    pImpl->peakValue = 0.0;
    pImpl->currentDrawdown = 0.0;
}

TradingStrategyExecutionModel::~TradingStrategyExecutionModel() = default;

bool TradingStrategyExecutionModel::validateSignal(const TradeSignal& signal) {
    if (signal.confidence < m_config.minConfidence) {
        return false;
    }

    if (signal.quantity <= 0 || signal.price <= 0) {
        return false;
    }

    if (pImpl->openPositions.size() >= m_config.maxOpenPositions &&
        signal.action == TradeSignal::Action::BUY) {
        return false;
    }

    return checkRiskLimits(signal);
}

double TradingStrategyExecutionModel::calculateOptimalPositionSize(
    const TradeSignal& signal) {
    
    // Kelly Criterion for position sizing
    double winProbability = signal.confidence;
    double winLossRatio = m_config.takeProfitPercent / m_config.stopLossPercent;
    
    double kellyFraction = winProbability - (1 - winProbability) / winLossRatio;
    kellyFraction = std::max(0.0, std::min(kellyFraction, 1.0));
    
    // Apply maximum position size constraint
    double maxSize = m_config.maxPositionSize;
    return std::min(maxSize, signal.quantity * kellyFraction);
}

bool TradingStrategyExecutionModel::checkRiskLimits(const TradeSignal& signal) {
    double potentialExposure = pImpl->totalExposure;
    
    if (signal.action == TradeSignal::Action::BUY) {
        potentialExposure += signal.quantity * signal.price;
    } else if (signal.action == TradeSignal::Action::SELL) {
        potentialExposure -= signal.quantity * signal.price;
    }

    // Check maximum drawdown
    double potentialDrawdown = (pImpl->peakValue - potentialExposure) 
                             / pImpl->peakValue;
    if (potentialDrawdown > m_config.maxDrawdown) {
        return false;
    }

    return true;
}

TradingStrategyExecutionModel::ExecutionResult 
TradingStrategyExecutionModel::executeStrategy(const TradeSignal& signal) {
    ExecutionResult result{};
    
    try {
        if (!validateSignal(signal)) {
            result.success = false;
            result.error = "Signal validation failed";
            return result;
        }

        double optimalSize = calculateOptimalPositionSize(signal);
        
        // Create Solana transaction
        SolanaFramework::Transaction tx{};
        tx.programId = "DEX_PROGRAM_ID"; // Replace with actual program ID
        
        if (signal.action == TradeSignal::Action::BUY) {
            tx.amount = static_cast<uint64_t>(optimalSize * signal.price);
            // Set up buy transaction details
        } else if (signal.action == TradeSignal::Action::SELL) {
            tx.amount = static_cast<uint64_t>(optimalSize * signal.price);
            // Set up sell transaction details
        }

        // Execute transaction
        if (m_solanaFramework->executeTransaction(tx)) {
            result.success = true;
            result.executedPrice = signal.price;
            result.executedQuantity = optimalSize;
            result.transactionId = tx.signature;
            
            // Update position tracking
            if (signal.action == TradeSignal::Action::BUY) {
                pImpl->totalExposure += optimalSize * signal.price;
            } else if (signal.action == TradeSignal::Action::SELL) {
                pImpl->totalExposure -= optimalSize * signal.price;
            }
            
            // Update peak value and drawdown
            pImpl->peakValue = std::max(pImpl->peakValue, pImpl->totalExposure);
            pImpl->currentDrawdown = (pImpl->peakValue - pImpl->totalExposure) 
                                   / pImpl->peakValue;
            
            pImpl->openPositions.push_back(result);
        } else {
            result.success = false;
            result.error = "Transaction execution failed";
        }

    } catch (const std::exception& e) {
        result.success = false;
        result.error = e.what();
    }

    return result;
}

void TradingStrategyExecutionModel::updateRiskParameters(
    const StrategyConfig& newConfig) {
    m_config = newConfig;
}

std::vector<TradingStrategyExecutionModel::ExecutionResult> 
TradingStrategyExecutionModel::getOpenPositions() {
    return pImpl->openPositions;
}

} // namespace xsol