// src/models/wallet_cluster_analysis_model.cpp

#include "wallet_cluster_analysis_model.h"
#include <algorithm>
#include <numeric>
#include <set>
#include <queue>

namespace xsol {


void WalletClusterAnalysisModel::buildTransactionGraph(
    const std::vector<SolanaFramework::Transaction>& transactions) {
    
    std::set<std::string> uniqueAddresses;
    
    // First pass: collect unique addresses
    for (const auto& tx : transactions) {
        uniqueAddresses.insert(tx.from);
        uniqueAddresses.insert(tx.to);
    }
    
    // Initialize adjacency matrix
    size_t n = uniqueAddresses.size();
    pImpl->transactionMatrix = graph::AdjacencyMatrix(n);
    
    // Create address to index mapping
    std::map<std::string, size_t> addressIndex;
    size_t idx = 0;
    for (const auto& addr : uniqueAddresses) {
        addressIndex[addr] = idx++;
    }
    
    // Build transaction matrix
    for (const auto& tx : transactions) {
        size_t fromIdx = addressIndex[tx.from];
        size_t toIdx = addressIndex[tx.to];
        pImpl->transactionMatrix.addEdge(fromIdx, toIdx, tx.amount);
        pImpl->totalVolume += tx.amount;
    }
}

double WalletClusterAnalysisModel::calculateWalletInfluence(
    const WalletNode& wallet) {
    
    double volumeWeight = 0.0;
    double connectivityWeight = 0.0;
    double balanceWeight = 0.0;
    
    // Calculate normalized weights
    if (pImpl->totalVolume > 0) {
        volumeWeight = wallet.activity / pImpl->totalVolume;
    }
    
    size_t maxConnections = pImpl->walletGraph.size() - 1;
    if (maxConnections > 0) {
        connectivityWeight = static_cast<double>(wallet.connections.size()) 
                           / maxConnections;
    }
    
    // Find maximum balance for normalization
    double maxBalance = 0.0;
    for (const auto& [_, node] : pImpl->walletGraph) {
        maxBalance = std::max(maxBalance, node.balance);
    }
    
    if (maxBalance > 0) {
        balanceWeight = wallet.balance / maxBalance;
    }
    
    // Combine metrics with weightings
    return 0.4 * volumeWeight + 
           0.3 * connectivityWeight + 
           0.3 * balanceWeight;
}

std::vector<std::vector<WalletClusterAnalysisModel::WalletNode>> 
WalletClusterAnalysisModel::detectCommunities(
    const std::vector<WalletNode>& wallets) {
    
    // Implement Louvain method for community detection
    std::vector<std::vector<WalletNode>> communities;
    std::vector<int> assignments(wallets.size(), -1);
    double modularity = 0.0;
    
    // Phase 1: Optimize modularity locally
    bool improved = true;
    while (improved) {
        improved = false;
        for (size_t i = 0; i < wallets.size(); ++i) {
            int bestCommunity = assignments[i];
            double bestModularity = modularity;
            
            // Try moving node i to each community
            std::set<int> neighborCommunities;
            for (const auto& neighbor : wallets[i].connections) {
                for (size_t j = 0; j < wallets.size(); ++j) {
                    if (wallets[j].address == neighbor) {
                        neighborCommunities.insert(assignments[j]);
                        break;
                    }
                }
            }
            
            for (int community : neighborCommunities) {
                assignments[i] = community;
                double newModularity = graph::calculateModularity(
                    pImpl->transactionMatrix, 
                    assignments
                );
                
                if (newModularity > bestModularity) {
                    bestModularity = newModularity;
                    bestCommunity = community;
                    improved = true;
                }
            }
            
            assignments[i] = bestCommunity;
            modularity = bestModularity;
        }
    }
    
    // Collect communities
    std::map<int, std::vector<WalletNode>> communityMap;
    for (size_t i = 0; i < wallets.size(); ++i) {
        communityMap[assignments[i]].push_back(wallets[i]);
    }
    
    for (const auto& [_, community] : communityMap) {
        communities.push_back(community);
    }
    
    return communities;
}

WalletClusterAnalysisModel::ClusterMetrics 
WalletClusterAnalysisModel::analyzeClusters(
    const std::vector<std::string>& addresses,
    const std::vector<SolanaFramework::Transaction>& transactions) {
    
    ClusterMetrics metrics{};
    pImpl->reset();
    
    // Build transaction graph
    buildTransactionGraph(transactions);
    
    // Create wallet nodes
    std::vector<WalletNode> wallets;
    for (const auto& addr : addresses) {
        WalletNode node;
        node.address = addr;
        
        // Analyze transactions for this wallet
        node.activity = 0.0;
        for (const auto& tx : transactions) {
            if (tx.from == addr || tx.to == addr) {
                node.activity += tx.amount;
                if (tx.from == addr) {
                    node.connections.push_back(tx.to);
                } else {
                    node.connections.push_back(tx.from);
                }
            }
        }
        
        // Remove duplicate connections
        std::sort(node.connections.begin(), node.connections.end());
        node.connections.erase(
            std::unique(node.connections.begin(), node.connections.end()),
            node.connections.end()
        );
        
        wallets.push_back(node);
        pImpl->walletGraph[addr] = node;
    }
    
    // Detect communities
    auto communities = detectCommunities(wallets);
    
    // Calculate metrics
    metrics.density = static_cast<double>(pImpl->transactionMatrix.getEdgeCount()) 
                     / (addresses.size() * (addresses.size() - 1));
    
    // Calculate centralization
    std::vector<double> influences;
    for (const auto& wallet : wallets) {
        double influence = calculateWalletInfluence(wallet);
        influences.push_back(influence);
        
        if (influence > m_config.significanceThreshold) {
            metrics.keyPlayers.push_back(wallet.address);
        }
    }
    
    // Calculate centralization using Gini coefficient
    metrics.centralization = graph::calculateGiniCoefficient(influences);
    
    // Calculate average balance and volume concentration
    double totalBalance = 0.0;
    for (const auto& wallet : wallets) {
        totalBalance += wallet.balance;
    }
    metrics.averageBalance = totalBalance / wallets.size();
    
    // Sort key players by influence
    std::sort(metrics.keyPlayers.begin(), metrics.keyPlayers.end(),
        [this](const std::string& a, const std::string& b) {
            return calculateWalletInfluence(pImpl->walletGraph[a]) >
                   calculateWalletInfluence(pImpl->walletGraph[b]);
        }
    );
    
    // Calculate volume concentration (Herfindahl-Hirschman Index)
    metrics.volumeConcentration = 0.0;
    for (const auto& wallet : wallets) {
        double marketShare = wallet.activity / pImpl->totalVolume;
        metrics.volumeConcentration += marketShare * marketShare;
    }
    
    return metrics;
}

std::vector<WalletClusterAnalysisModel::WalletNode>
WalletClusterAnalysisModel::findInfluentialWallets(
    const std::vector<std::string>& addresses) {
    
    std::vector<WalletNode> influential;
    
    for (const auto& addr : addresses) {
        auto it = pImpl->walletGraph.find(addr);
        if (it != pImpl->walletGraph.end() && 
            calculateWalletInfluence(it->second) > m_config.significanceThreshold) {
            influential.push_back(it->second);
        }
    }
    
    // Sort by influence
    std::sort(influential.begin(), influential.end(),
        [this](const WalletNode& a, const WalletNode& b) {
            return calculateWalletInfluence(a) > calculateWalletInfluence(b);
        }
    );
    
    return influential;
}

void WalletClusterAnalysisModel::updateConfig(const AnalysisConfig& newConfig) {
    m_config = newConfig;
}

} // namespace xsol
