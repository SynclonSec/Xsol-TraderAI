#include "solana_framework.h"
#include <iostream>
#include <curl/curl.h>
#include <nlohmann/json.hpp>

namespace xsol
{

    struct SolanaFramework::Impl
    {
        CURL *curl;
        std::string api_key;
        std::string endpoint;
        std::string currency;
    };

    SolanaFramework::SolanaFramework() : pImpl(std::make_unique<Impl>())
    {
        // fetch Solana API info.
        const char *api_key = std::getenv("SOLANA_API_KEY");
        const char *endpoint = std::getenv("SOLANA_ENDPOINT");
        const char *currency = std::getenv("SOLANA_CURRENCY");
        if (!api_key || !endpoint || !currency)
        {
            throw std::runtime_error("Environment variables SOLANA_API_KEY and SOLANA_ENDPOINT must be set");
        }
        pImpl->api_key = api_key;
        pImpl->endpoint = endpoint;
        pImpl->currency = currency;

        // Initialize CURL
        pImpl->curl = curl_easy_init();
        if (!pImpl->curl)
        {
            throw std::runtime_error("Failed to initialize CURL");
        }
    }

    SolanaFramework::~SolanaFramework()
    {
        if (pImpl->curl)
        {
            curl_easy_cleanup(pImpl->curl);
        }
    }

    bool SolanaFramework::initializeSolanaConnection(const std::string &)
    {
        // todo: websocket connection for events
        return true;
    }

    bool SolanaFramework::executeTransaction(const Transaction &tx)
    {
        // Execute the transaction using the Solana SDK
        // Example: return solana::rpc::sendTransaction(tx.signature, tx.from, tx.to, tx.amount, tx.programId);
        std::cout << "Executing transaction from " << tx.from << " to " << tx.to << " amount " << tx.amount << std::endl;
        return true;
    }

    SolanaFramework::MarketData SolanaFramework::getMarketData(const std::string &market)
    {
        MarketData data;
        // Fetch market data using the Solana SDK
        // Example: data = solana::rpc::getMarketData(market);
        std::cout << "Fetching market data for " << market << std::endl;
        return data;
    }

    double SolanaFramework::getAccountBalance(const std::string &address)
    {
        std::cout << "Fetching account balance for " << address << std::endl;

        // Construct API request
        std::string url = pImpl->endpoint + address + "/balances?currencies=" + pImpl->currency;
        struct curl_slist *headers = nullptr;
        std::string api_key_header = "X-API-KEY: " + pImpl->api_key;
        headers = curl_slist_append(headers, api_key_header.c_str());

        curl_easy_setopt(pImpl->curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(pImpl->curl, CURLOPT_HTTPHEADER, headers);
        curl_easy_setopt(pImpl->curl, CURLOPT_WRITEFUNCTION, [](void *contents, size_t size, size_t nmemb, std::string *s) -> size_t
                         {
        size_t newLength = size * nmemb;
        try {
            s->append((char*)contents, newLength);
        } catch (std::bad_alloc& e) {
            // Handle memory problem
            return 0;
        }
        return newLength; });

        // Parse balance.
        std::string response;
        curl_easy_setopt(pImpl->curl, CURLOPT_WRITEDATA, &response);
        CURLcode res = curl_easy_perform(pImpl->curl);
        if (res != CURLE_OK)
        {
            std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << std::endl;
            return 0.0;
        }
        curl_slist_free_all(headers);
        auto json_response = nlohmann::json::parse(response);
        double balance = 0.0;
        for (const auto &item : json_response)
        {
            if (item["currency"] == "usdc")
            {
                balance = std::stod(item["balances"]["total"].get<std::string>());
                break;
            }
        }

        return balance;
    }

    std::vector<SolanaFramework::Transaction> SolanaFramework::getRecentTransactions(const std::string &address, uint32_t limit)
    {
        std::vector<Transaction> transactions;
        // Fetch recent transactions using the Solana SDK
        // Example: transactions = solana::rpc::getRecentTransactions(address, limit);
        std::cout << "Fetching recent transactions for " << address << " with limit " << limit << std::endl;
        return transactions;
    }

} // namespace xsol