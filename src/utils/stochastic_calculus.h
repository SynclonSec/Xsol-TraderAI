#ifndef STOCHASTIC_CALCULUS_H
#define STOCHASTIC_CALCULUS_H

#include <vector>
#include <functional>

typedef std::vector<double> (*SDEFunction)(double, const std::vector<double>&);

class StochasticCalculus {
public:
    StochasticCalculus(double dt);

    std::vector<double> generateBrownianMotion(double t0, double tf, double initialValue);
    std::vector<double> simulateStochasticProcess(double t0, double tf, const std::vector<double>& y0,
                                                  SDEFunction drift, SDEFunction diffusion);
    double calculateExpectation(std::function<double(const std::vector<double>&)> func,
                                const std::vector<std::vector<double>>& samples);
    std::vector<std::vector<double>> generateSamples(double t0, double tf, const std::vector<double>& y0,
                                                     SDEFunction drift, SDEFunction diffusion, int numSamples);
    double calculateVariance(std::function<double(const std::vector<double>&)> func,
                              const std::vector<std::vector<double>>& samples);
    double calculateCovariance(std::function<double(const std::vector<double>&)> func1,
                                std::function<double(const std::vector<double>&)> func2,
                                const std::vector<std::vector<double>>& samples);
    double calculateCorrelation(std::function<double(const std::vector<double>&)> func1,
                                std::function<double(const std::vector<double>&)> func2,
                                const std::vector<std::vector<double>>& samples);
    std::vector<double> generateGeometricBrownianMotion(double t0, double tf, double initialValue, double mu, double sigma);
    std::vector<double> generateOrnsteinUhlenbeckProcess(double t0, double tf, double initialValue, double theta, double mu, double sigma);
    std::vector<int> generatePoissonProcess(double t0, double tf, double lambda);

private:
    double dt;
    std::mt19937 generator;

    std::vector<double> addVectors(const std::vector<double>& a, const std::vector<double>& b);
    std::vector<double> scaleVector(const std::vector<double>& v, double scalar);
    double generateBrownianIncrement();
    std::vector<double> generateBrownianIncrementVector(size_t size);
};

std::vector<double> driftFunction(double t, const std::vector<double>& y);
std::vector<double> diffusionFunction(double t, const std::vector<double>& y);
double exampleFunction(const std::vector<double>& y);

#endif // STOCHASTIC_CALCULUS_H
