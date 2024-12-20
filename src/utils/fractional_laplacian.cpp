#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <complex>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <omp.h>

void log(const std::string& message) {
    std::ofstream logFile("fractional_laplacian.log", std::ios_base::app);
    if (logFile.is_open()) {
        logFile << message << std::endl;
        logFile.close();
    }
}

std::vector<double> fractionalLaplacian(const std::vector<double>& u, double alpha, double h) {
    int N = u.size();   
    std::vector<double> result(N, 0.0);

    if (alpha <= 0 || alpha >= 2) {
        log("Error: alpha must be in the range (0, 2).");
        throw std::invalid_argument("alpha must be in the range (0, 2).");
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                double diff = std::abs(i - j) * h;
                sum += (u[i] - u[j]) / std::pow(diff, 1 + alpha);
            }
        }
        result[i] = sum * -std::tgamma(2 - alpha) / (2 * std::tgamma((1 - alpha) / 2) * std::tgamma((1 + alpha) / 2));
    }

    return result;
}

void displayVector(const std::vector<double>& vec) {
    for (const auto& value : vec) {
        std::cout << std::fixed << std::setprecision(4) << value << " ";
    }
    std::cout << std::endl;
}

double l2Norm(const std::vector<double>& vec) {
    double sum = 0.0;
    for (const auto& value : vec) {
        sum += value * value;
    }
    return std::sqrt(sum);
}

double mean(const std::vector<double>& vec) {
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    return sum / vec.size();
}

double standardDeviation(const std::vector<double>& vec) {
    double vecMean = mean(vec);
    double sum = 0.0;
    for (const auto& value : vec) {
        sum += (value - vecMean) * (value - vecMean);
    }
    return std::sqrt(sum / vec.size());
}

std::vector<double> simpleMovingAverage(const std::vector<double>& vec, int windowSize) {
    if (windowSize <= 0 || windowSize > vec.size()) {
        log("Error: Invalid window size for simple moving average.");
        throw std::invalid_argument("Invalid window size for simple moving average.");
    }

    std::vector<double> result(vec.size(), 0.0);
    for (size_t i = 0; i < vec.size(); ++i) {
        int start = std::max(0, static_cast<int>(i) - windowSize + 1);
        int end = i + 1;
        double sum = std::accumulate(vec.begin() + start, vec.begin() + end, 0.0);
        result[i] = sum / (end - start);
    }
    return result;
}

std::vector<double> exponentialMovingAverage(const std::vector<double>& vec, double smoothingFactor) {
    if (smoothingFactor <= 0 || smoothingFactor >= 1) {
        log("Error: Smoothing factor must be in the range (0, 1).");
        throw std::invalid_argument("Smoothing factor must be in the range (0, 1).");
    }

    std::vector<double> result(vec.size(), 0.0);
    result[0] = vec[0];
    for (size_t i = 1; i < vec.size(); ++i) {
        result[i] = smoothingFactor * vec[i] + (1 - smoothingFactor) * result[i - 1];
    }
    return result;
}

std::vector<double> fractionalDerivative(const std::vector<double>& u, double alpha, double h) {
    int N = u.size();
    std::vector<double> result(N, 0.0);

    if (alpha <= 0 || alpha >= 1) {
        log("Error: alpha must be in the range (0, 1) for fractional derivatives.");
        throw std::invalid_argument("alpha must be in the range (0, 1) for fractional derivatives.");
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                double diff = std::abs(i - j) * h;
                sum += (u[i] - u[j]) / std::pow(diff, alpha);
            }
        }
        result[i] = sum * std::tgamma(1 - alpha);
    }

    return result;
}

std::vector<double> fractionalIntegral(const std::vector<double>& u, double alpha, double h) {
    int N = u.size();
    std::vector<double> result(N, 0.0);

    if (alpha <= 0 || alpha >= 1) {
        log("Error: alpha must be in the range (0, 1) for fractional integrals.");
        throw std::invalid_argument("alpha must be in the range (0, 1) for fractional integrals.");
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                double diff = std::abs(i - j) * h;
                sum += (u[i] - u[j]) * std::pow(diff, alpha - 1);
            }
        }
        result[i] = sum * std::tgamma(alpha);
    }

    return result;
}

std::vector<double> fractionalDiffusion(const std::vector<double>& u, double alpha, double h, double dt) {
    int N = u.size();
    std::vector<double> result(N, 0.0);

    if (alpha <= 0 || alpha >= 2) {
        log("Error: alpha must be in the range (0, 2) for fractional diffusion.");
        throw std::invalid_argument("alpha must be in the range (0, 2) for fractional diffusion.");
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                double diff = std::abs(i - j) * h;
                sum += (u[i] - u[j]) / std::pow(diff, 1 + alpha);
            }
        }
        result[i] = sum * -std::tgamma(2 - alpha) / (2 * std::tgamma((1 - alpha) / 2) * std::tgamma((1 + alpha) / 2)) * dt;
    }

    return result;
}

std::vector<double> fractionalBrownianMotion(int N, double H) {
    std::vector<double> result(N, 0.0);

    if (H <= 0 || H >= 1) {
        log("Error: H must be in the range (0, 1) for fractional Brownian motion.");
        throw std::invalid_argument("H must be in the range (0, 1) for fractional Brownian motion.");
    }

    for (int i = 1; i < N; ++i) {
        result[i] = result[i - 1] + std::pow(i, H - 0.5) * (std::rand() / static_cast<double>(RAND_MAX) - 0.5);
    }

    return result;
}

void runUnitTests() {
    try {
        std::vector<double> u = {1, 2, 3, 4, 5};
        double alpha = 1.5;
        double h = 1.0;
        std::vector<double> result = fractionalLaplacian(u, alpha, h);

        std::vector<double> expected = {-0.3333, -0.1089, 0.0000, 0.1089, 0.3333};
        for (size_t i = 0; i < result.size(); ++i) {
            if (std::abs(result[i] - expected[i]) >= 1e-4) {
                log("Unit test failed: result[" + std::to_string(i) + "] = " + std::to_string(result[i]) +
                    ", expected[" + std::to_string(i) + "] = " + std::to_string(expected[i]));
                assert(std::abs(result[i] - expected[i]) < 1e-4);
            }
        }

        std::cout << "All unit tests passed." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Unit test failed: " << e.what() << std::endl;
    }
}

template <typename Func, typename... Args>
std::vector<double> measureExecutionTime(Func func, Args&&... args) {
    auto start = std::chrono::high_resolution_clock::now();
    auto result = func(std::forward<Args>(args)...);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;
    return result;
}

int main() {
    try {
        std::vector<double> u = {1, 2, 3, 4, 5};
        double alpha = 1.5;
        double h = 1.0;

        std::vector<double> result = measureExecutionTime(fractionalLaplacian, u, alpha, h);

        std::cout << "Fractional Laplacian Result:" << std::endl;
        displayVector(result);

        std::cout << "L2 Norm: " << l2Norm(result) << std::endl;
        std::cout << "Mean: " << mean(result) << std::endl;
        std::cout << "Standard Deviation: " << standardDeviation(result) << std::endl;

        int windowSize = 3;
        std::vector<double> smaResult = simpleMovingAverage(result, windowSize);
        std::cout << "Simple Moving Average (window size = " << windowSize << "):" << std::endl;
        displayVector(smaResult);

        double smoothingFactor = 0.5;
        std::vector<double> emaResult = exponentialMovingAverage(result, smoothingFactor);
        std::cout << "Exponential Moving Average (smoothing factor = " << smoothingFactor << "):" << std::endl;
        displayVector(emaResult);

        std::vector<double> fdResult = fractionalDerivative(result, 0.5, h);
        std::cout << "Fractional Derivative Result:" << std::endl;
        displayVector(fdResult);

        std::vector<double> fiResult = fractionalIntegral(result, 0.5, h);
        std::cout << "Fractional Integral Result:" << std::endl;
        displayVector(fiResult);

        std::vector<double> fdifResult = fractionalDiffusion(result, alpha, h, 0.1);
        std::cout << "Fractional Diffusion Result:" << std::endl;
        displayVector(fdifResult);

        std::vector<double> fbmResult = fractionalBrownianMotion(10, 0.5);
        std::cout << "Fractional Brownian Motion Result:" << std::endl;
        displayVector(fbmResult);

        runUnitTests();
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Unexpected error: " << e.what() << std::endl;
    }

    return 0;
}
