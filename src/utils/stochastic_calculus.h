
#pragma once

#include <vector>
#include <random>
#include <cmath>

namespace xsol {
namespace math {

class StochasticCalculus {
public:
    struct BrownianMotion {
        double drift;
        double volatility;
        std::vector<double> path;
        double timeStep;
    };

    struct HestonParameters {
        double kappa;    // Mean reversion speed
        double theta;    // Long-run variance
        double sigma;    // Volatility of variance
        double rho;      // Correlation
        double v0;       // Initial variance
    };

    StochasticCalculus(unsigned seed = std::random_device{}());

    // Generate Brownian motion path
    BrownianMotion generateBrownianMotion(
        double T,           // Time horizon
        size_t steps,       // Number of steps
        double drift = 0.0, // Drift term
        double vol = 1.0    // Volatility
    );

    // Implement Heston stochastic volatility model
    std::vector<double> simulateHestonModel(
        const HestonParameters& params,
        double S0,    // Initial price
        double T,     // Time horizon
        size_t steps  // Number of steps
    );

    // Calculate Ito integral
    double calculateItoIntegral(
        const std::vector<double>& process,
        const std::vector<double>& integrand,
        double dt
    );

private:
    std::mt19937_64 m_rng;
    std::normal_distribution<double> m_normalDist;
};

} // namespace math
} // namespace xsol