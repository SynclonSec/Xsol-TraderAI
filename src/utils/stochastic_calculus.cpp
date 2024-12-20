#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <functional>

typedef std::vector<double> (*SDEFunction)(double, const std::vector<double>&);

class StochasticCalculus {
public:
    StochasticCalculus(double dt) : dt(dt) {
        std::random_device rd;
        generator.seed(rd());
    }

    std::vector<double> generateBrownianMotion(double t0, double tf, double initialValue) {
        int n = static_cast<int>((tf - t0) / dt);
        std::vector<double> brownianMotion(n + 1, initialValue);
        for (int i = 1; i <= n; ++i) {
            brownianMotion[i] = brownianMotion[i - 1] + generateBrownianIncrement();
        }
        return brownianMotion;
    }

    std::vector<double> simulateStochasticProcess(double t0, double tf, const std::vector<double>& y0,
                                                  SDEFunction drift, SDEFunction diffusion) {
        int n = static_cast<int>((tf - t0) / dt);
        std::vector<double> y = y0;
        double t = t0;

        for (int i = 0; i < n; ++i) {
            std::vector<double> driftTerm = drift(t, y);
            std::vector<double> diffusionTerm = diffusion(t, y);
            std::vector<double> dW = generateBrownianIncrementVector(y.size());

            y = addVectors(y, scaleVector(driftTerm, dt));
            y = addVectors(y, scaleVector(diffusionTerm, std::sqrt(dt) * dW[0]));

            t += dt;
        }

        return y;
    }

    double calculateExpectation(std::function<double(const std::vector<double>&)> func,
                                const std::vector<std::vector<double>>& samples) {
        double sum = 0.0;
        for (const auto& sample : samples) {
            sum += func(sample);
        }
        return sum / samples.size();
    }

    std::vector<std::vector<double>> generateSamples(double t0, double tf, const std::vector<double>& y0,
                                                     SDEFunction drift, SDEFunction diffusion, int numSamples) {
        std::vector<std::vector<double>> samples(numSamples);
        for (int i = 0; i < numSamples; ++i) {
            samples[i] = simulateStochasticProcess(t0, tf, y0, drift, diffusion);
        }
        return samples;
    }

    double calculateVariance(std::function<double(const std::vector<double>&)> func,
                              const std::vector<std::vector<double>>& samples) {
        double mean = calculateExpectation(func, samples);
        double sumSquaredDiff = 0.0;
        for (const auto& sample : samples) {
            double diff = func(sample) - mean;
            sumSquaredDiff += diff * diff;
        }
        return sumSquaredDiff / samples.size();
    }

    double calculateCovariance(std::function<double(const std::vector<double>&)> func1,
                                std::function<double(const std::vector<double>&)> func2,
                                const std::vector<std::vector<double>>& samples) {
        double mean1 = calculateExpectation(func1, samples);
        double mean2 = calculateExpectation(func2, samples);
        double sumProductDiff = 0.0;
        for (const auto& sample : samples) {
            double diff1 = func1(sample) - mean1;
            double diff2 = func2(sample) - mean2;
            sumProductDiff += diff1 * diff2;
        }
        return sumProductDiff / samples.size();
    }

    double calculateCorrelation(std::function<double(const std::vector<double>&)> func1,
                                std::function<double(const std::vector<double>&)> func2,
                                const std::vector<std::vector<double>>& samples) {
        double cov = calculateCovariance(func1, func2, samples);
        double var1 = calculateVariance(func1, samples);
        double var2 = calculateVariance(func2, samples);
        return cov / std::sqrt(var1 * var2);
    }

    std::vector<double> generateGeometricBrownianMotion(double t0, double tf, double initialValue, double mu, double sigma) {
        int n = static_cast<int>((tf - t0) / dt);
        std::vector<double> gbm(n + 1, initialValue);
        for (int i = 1; i <= n; ++i) {
            double dW = generateBrownianIncrement();
            gbm[i] = gbm[i - 1] * std::exp((mu - 0.5 * sigma * sigma) * dt + sigma * dW);
        }
        return gbm;
    }

    std::vector<double> generateOrnsteinUhlenbeckProcess(double t0, double tf, double initialValue, double theta, double mu, double sigma) {
        int n = static_cast<int>((tf - t0) / dt);
        std::vector<double> ou(n + 1, initialValue);
        for (int i = 1; i <= n; ++i) {
            double dW = generateBrownianIncrement();
            ou[i] = ou[i - 1] + theta * (mu - ou[i - 1]) * dt + sigma * dW;
        }
        return ou;
    }

    std::vector<int> generatePoissonProcess(double t0, double tf, double lambda) {
        int n = static_cast<int>((tf - t0) / dt);
        std::vector<int> poissonProcess(n + 1, 0);
        std::poisson_distribution<int> distribution(lambda * dt);
        for (int i = 1; i <= n; ++i) {
            poissonProcess[i] = poissonProcess[i - 1] + distribution(generator);
        }
        return poissonProcess;
    }

private:
    double dt;
    std::mt19937 generator;

    std::vector<double> addVectors(const std::vector<double>& a, const std::vector<double>& b) {
        std::vector<double> result(a.size());
        for (size_t i = 0; i < a.size(); ++i) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    std::vector<double> scaleVector(const std::vector<double>& v, double scalar) {
        std::vector<double> result(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            result[i] = v[i] * scalar;
        }
        return result;
    }

    double generateBrownianIncrement() {
        std::normal_distribution<double> distribution(0.0, std::sqrt(dt));
        return distribution(generator);
    }

    std::vector<double> generateBrownianIncrementVector(size_t size) {
        std::vector<double> dW(size);
        for (size_t i = 0; i < size; ++i) {
            dW[i] = generateBrownianIncrement();
        }
        return dW;
    }
};

std::vector<double> driftFunction(double t, const std::vector<double>& y) {
    double mu = 0.1;
    return {mu * y[0]};
}

std::vector<double> diffusionFunction(double t, const std::vector<double>& y) {
    double sigma = 0.2;
    return {sigma * y[0]};
}

double exampleFunction(const std::vector<double>& y) {
    return y[0];
}

int main() {
    double t0 = 0.0;
    double tf = 1.0;
    std::vector<double> y0 = {1.0};
    double dt = 0.01;
    int numSamples = 1000;

    StochasticCalculus sc(dt);

    std::vector<double> brownianMotion = sc.generateBrownianMotion(t0, tf, 0.0);
    std::cout << "Brownian Motion at t = " << tf << ": " << brownianMotion.back() << std::endl;

    std::vector<double> result = sc.simulateStochasticProcess(t0, tf, y0, driftFunction, diffusionFunction);
    std::cout << "Stochastic Process at t = " << tf << ": y = " << result[0] << std::endl;

    std::vector<std::vector<double>> samples = sc.generateSamples(t0, tf, y0, driftFunction, diffusionFunction, numSamples);

    double expectation = sc.calculateExpectation(exampleFunction, samples);
    std::cout << "Expectation of the process at t = " << tf << ": E[Y] = " << expectation << std::endl;

    double variance = sc.calculateVariance(exampleFunction, samples);
    std::cout << "Variance of the process at t = " << tf << ": Var[Y] = " << variance << std::endl;

    std::vector<double> gbm = sc.generateGeometricBrownianMotion(t0, tf, 1.0, 0.1, 0.2);
    std::cout << "Geometric Brownian Motion at t = " << tf << ": " << gbm.back() << std::endl;

    std::vector<double> ou = sc.generateOrnsteinUhlenbeckProcess(t0, tf, 0.0, 0.5, 1.0, 0.2);
    std::cout << "Ornstein-Uhlenbeck Process at t = " << tf << ": " << ou.back() << std::endl;

    std::vector<int> poissonProcess = sc.generatePoissonProcess(t0, tf, 1.0);
    std::cout << "Poisson Process at t = " << tf << ": " << poissonProcess.back() << std::endl;

    return 0;
}
