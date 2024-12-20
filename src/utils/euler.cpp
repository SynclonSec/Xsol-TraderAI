#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <cassert>

void log(const std::string& message) {
    std::ofstream logFile("trading_agent.log", std::ios_base::app);
    if (logFile.is_open()) {
        logFile << message << std::endl;
        logFile.close();
    }
}

std::complex<double> eulerFormula(double x) {
    if (std::isnan(x) || std::isinf(x)) {
        log("Error: Input value is not a valid real number.");
        throw std::invalid_argument("Input value is not a valid real number.");
    }
    std::complex<double> result = std::exp(std::complex<double>(0, x));
    log("Euler's formula result for x=" + std::to_string(x) + ": " + std::to_string(result.real()) + " + " + std::to_string(result.imag()) + "i");
    return result;
}

std::vector<std::complex<double>> dft(const std::vector<double>& priceData) {
    int N = priceData.size();
    std::vector<std::complex<double>> result(N);
    for (int k = 0; k < N; ++k) {
        result[k] = std::complex<double>(0, 0);
        for (int n = 0; n < N; ++n) {
            double angle = 2 * M_PI * k * n / N;
            result[k] += priceData[n] * eulerFormula(-angle);
        }
    }
    return result;
}

void displayComplexNumber(const std::complex<double>& number) {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Real part: " << number.real() << std::endl;
    std::cout << "Imaginary part: " << number.imag() << "i" << std::endl;
}

void runUnitTests() {
    try {
        std::complex<double> result1 = eulerFormula(0);
        assert(std::abs(result1.real() - 1.0) < 1e-9 && std::abs(result1.imag() - 0.0) < 1e-9);

        std::complex<double> result2 = eulerFormula(M_PI);
        assert(std::abs(result2.real() + 1.0) < 1e-9 && std::abs(result2.imag() - 0.0) < 1e-9);

        std::complex<double> result3 = eulerFormula(M_PI / 2);
        assert(std::abs(result3.real() - 0.0) < 1e-9 && std::abs(result3.imag() - 1.0) < 1e-9);

        std::cout << "All unit tests passed." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Unit test failed: " << e.what() << std::endl;
    }
}

int main() {
    try {
        std::vector<double> priceData = {100, 101, 102, 103, 104, 105, 106, 107, 108, 109};

        std::vector<std::complex<double>> dftResult = dft(priceData);

        std::cout << "DFT Result:" << std::endl;
        for (const auto& value : dftResult) {
            displayComplexNumber(value);
        }

        runUnitTests();
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Unexpected error: " << e.what() << std::endl;
    }

    return 0;
}

