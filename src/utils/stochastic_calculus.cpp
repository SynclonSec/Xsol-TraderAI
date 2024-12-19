
#include "stochastic_calculus.h"
#include <stdexcept>

namespace xsol {
namespace math {

StochasticCalculus::StochasticCalculus(unsigned seed)
    : m_rng(seed)
    , m_normalDist(0.0, 1.0) {}

StochasticCalculus::BrownianMotion 
StochasticCalculus::generateBrownianMotion(
    double T,
    size_t steps,
    double drift,
    double vol) {

    if (T <= 0 || steps == 0) {
        throw std::invalid_argument("Invalid parameters for Brownian motion");
    }

    double dt = T / steps;
    double sqrtDt = std::sqrt(dt);

    BrownianMotion bm{
        .drift = drift,
        .volatility = vol,
        .path = std::vector<double>(steps + 1, 0.0),
        .timeStep = dt
    };

    for (size_t i = 1; i <= steps; ++i) {
        double dW = m_normalDist(m_rng) * sqrtDt;
        bm.path[i] = bm.path[i-1] + drift * dt + vol * dW;
    }

    return bm;
}

std::vector<double> 
StochasticCalculus::simulateHestonModel(
    const HestonParameters& params,
    double S0,
    double T,
    size_t steps) {

    if (T <= 0 || steps == 0 || S0 <= 0) {
        throw std::invalid_argument("Invalid parameters for Heston model");
    }

    double dt = T / steps;
    double sqrtDt = std::sqrt(dt);

    std::vector<double> prices(steps + 1, S0);
    std::vector<double> variances(steps + 1, params.v0);

    for (size_t i = 1; i <= steps; ++i) {
        double z1 = m_normalDist(m_rng);
        double z2 = params.rho * z1 + 
                   std::sqrt(1 - params.rho * params.rho) * m_normalDist(m_rng);

        // Update variance process
        variances[i] = std::max(0.0, variances[i-1] + 
            params.kappa * (params.theta - variances[i-1]) * dt +
            params.sigma * std::sqrt(variances[i-1]) * z1 * sqrtDt);

        // Update price process
        prices[i] = prices[i-1] * std::exp(
            (- 0.5 * variances[i-1]) * dt +
            std::sqrt(variances[i-1]) * z2 * sqrtDt
        );
    }

    return prices;
}

double StochasticCalculus::calculateItoIntegral(
    const std::vector<double>& process,
    const std::vector<double>& integrand,
    double dt) {
    
    if (process.size() != integrand.size()) {
        throw std::invalid_argument("Process and integrand must have same size");
    }

    double integral = 0.0;
    for (size_t i = 0; i < process.size() - 1; ++i) {
        integral += integrand[i] * (process[i+1] - process[i]);
    }

    return integral;
}

} // namespace math
} // namespace xsol