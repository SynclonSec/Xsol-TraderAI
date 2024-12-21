
#pragma once

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include "mathematical_tools.h"

namespace xsol {

class SolanaFramework {
public:
    struct MarketData {
        double price;
        double volume;
        double volatility;
        int64_t timestamp;
        std::vector<double> orderBookBids;
        std::vector<double> orderBookAsks;
    };

    struct Transaction {
        std::string signature;
        std::string from;
        std::string to;
        uint64_t amount;
        std::string programId;
    };

    struct Order {
        enum class Action { BUY, SELL };
        Action action;
        double quantity;
        double price;
        double acceptableSlippage;
        std::string market;
    };

    SolanaFramework();
    ~SolanaFramework();

    bool initializeSolanaConnection(const std::string& endpoint);
    MarketData getMarketData(const std::string& market);
    double getAccountBalance(const std::string& address);
    std::vector<Transaction> getRecentTransactions(const std::string& address, uint32_t limit = 10);

    bool executeOrder(const Order& order);
private:
    struct Impl;
    std::unique_ptr<Impl> pImpl; // PIMPL idiom for ABI stability
};

} // namespace xsol
