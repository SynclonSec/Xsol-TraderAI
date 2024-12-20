
#pragma once

#include <vector>
#include <string>
#include <map>
#include <memory>
#include "../solana/solana_framework.h"

namespace xsol {

class WalletClusterAnalysisModel {
public:
    struct WalletNode {
        std::string address;
        double balance;
        std::vector<std::string> connections;
        double influence;
        double activity;
        std::map<std::string, double> tokenBalances;
    };

    struct ClusterMetrics {
        double density;
        double centralization;
        double averageBalance;
        double volumeConcentration;
        std::vector<std::string> keyPlayers;
    };

    struct AnalysisConfig {
        double minBalance;
        int minTransactions;
        double timeWindowHours;
        int maxClusters;
        double significanceThreshold;
    };

    WalletClusterAnalysisModel(const AnalysisConfig& config);
    ~WalletClusterAnalysisModel();

    // Analyze wallet clusters and return metrics
    ClusterMetrics analyzeClusters(
        const std::vector<std::string>& addresses,
        const std::vector<SolanaFramework::Transaction>& transactions
    );

    // Identify influential wallets
    std::vector<WalletNode> findInfluentialWallets(
        const std::vector<std::string>& addresses
    );

    // Update analysis parameters
    void updateConfig(const AnalysisConfig& newConfig);

private:
    struct Impl;
    std::unique_ptr<Impl> pImpl;
    AnalysisConfig m_config;

    // Helper methods
    double calculateWalletInfluence(const WalletNode& wallet);
    std::vector<std::vector<WalletNode>> detectCommunities(
        const std::vector<WalletNode>& wallets
    );
    void buildTransactionGraph(
        const std::vector<SolanaFramework::Transaction>& transactions
    );
};

} // namespace xsol