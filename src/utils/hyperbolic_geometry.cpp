#ifndef HYPERBOLIC_GEOMETRY_H
#define HYPERBOLIC_GEOMETRY_H

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <complex>
#include <iomanip>
#include <cassert>

// Function declarations
double hyperbolicSine(double x);
double hyperbolicCosine(double x);
double hyperbolicTangent(double x);
double hyperbolicCotangent(double x);
double hyperbolicSecant(double x);
double hyperbolicCosecant(double x);

std::complex<double> hyperbolicSineComplex(const std::complex<double>& z);
std::complex<double> hyperbolicCosineComplex(const std::complex<double>& z);
std::complex<double> hyperbolicTangentComplex(const std::complex<double>& z);
std::complex<double> hyperbolicCotangentComplex(const std::complex<double>& z);
std::complex<double> hyperbolicSecantComplex(const std::complex<double>& z);
std::complex<double> hyperbolicCosecantComplex(const std::complex<double>& z);

void displayComplexNumber(const std::complex<double>& number);
void runUnitTests();

// Hyperbolic Sine
double hyperbolicSine(double x) {
    return std::sinh(x);
}

// Hyperbolic Cosine
double hyperbolicCosine(double x) {
    return std::cosh(x);
}

// Hyperbolic Tangent
double hyperbolicTangent(double x) {
    return std::tanh(x);
}

// Hyperbolic Cotangent
double hyperbolicCotangent(double x) {
    if (std::abs(x) < 1e-9) {
        throw std::invalid_argument("Hyperbolic cotangent is undefined for x = 0.");
    }
    return 1.0 / std::tanh(x);
}

// Hyperbolic Secant
double hyperbolicSecant(double x) {
    return 1.0 / std::cosh(x);
}

// Hyperbolic Cosecant
double hyperbolicCosecant(double x) {
    if (std::abs(x) < 1e-9) {
        throw std::invalid_argument("Hyperbolic cosecant is undefined for x = 0.");
    }
    return 1.0 / std::sinh(x);
}

// Hyperbolic Sine for Complex Numbers
std::complex<double> hyperbolicSineComplex(const std::complex<double>& z) {
    return std::sinh(z);
}

// Hyperbolic Cosine for Complex Numbers
std::complex<double> hyperbolicCosineComplex(const std::complex<double>& z) {
    return std::cosh(z);
}

// Hyperbolic Tangent for Complex Numbers
std::complex<double> hyperbolicTangentComplex(const std::complex<double>& z) {
    return std::tanh(z);
}

// Hyperbolic Cotangent for Complex Numbers
std::complex<double> hyperbolicCotangentComplex(const std::complex<double>& z) {
    if (std::abs(z) < 1e-9) {
        throw std::invalid_argument("Hyperbolic cotangent is undefined for z = 0.");
    }
    return 1.0 / std::tanh(z);
}

// Hyperbolic Secant for Complex Numbers
std::complex<double> hyperbolicSecantComplex(const std::complex<double>& z) {
    return 1.0 / std::cosh(z);
}

// Hyperbolic Cosecant for Complex Numbers
std::complex<double> hyperbolicCosecantComplex(const std::complex<double>& z) {
    if (std::abs(z) < 1e-9) {
        throw std::invalid_argument("Hyperbolic cosecant is undefined for z = 0.");
    }
    return 1.0 / std::sinh(z);
}

// Display Complex Number
void displayComplexNumber(const std::complex<double>& number) {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Real part: " << number.real() << std::endl;
    std::cout << "Imaginary part: " << number.imag() << "i" << std::endl;
}

// Run Unit Tests
void runUnitTests() {
    try {
        // Test hyperbolicSine
        double result1 = hyperbolicSine(0);
        assert(std::abs(result1 - 0.0) < 1e-9);

        double result2 = hyperbolicSine(1);
        assert(std::abs(result2 - std::sinh(1)) < 1e-9);

        // Test hyperbolicCosine
        double result3 = hyperbolicCosine(0);
        assert(std::abs(result3 - 1.0) < 1e-9);

        double result4 = hyperbolicCosine(1);
        assert(std::abs(result4 - std::cosh(1)) < 1e-9);

        // Test hyperbolicTangent
        double result5 = hyperbolicTangent(0);
        assert(std::abs(result5 - 0.0) < 1e-9);

        double result6 = hyperbolicTangent(1);
        assert(std::abs(result6 - std::tanh(1)) < 1e-9);

        // Test hyperbolicCotangent
        double result7 = hyperbolicCotangent(1);
        assert(std::abs(result7 - (1.0 / std::tanh(1))) < 1e-9);

        // Test hyperbolicSecant
        double result8 = hyperbolicSecant(0);
        assert(std::abs(result8 - 1.0) < 1e-9);

        double result9 = hyperbolicSecant(1);
        assert(std::abs(result9 - (1.0 / std::cosh(1))) < 1e-9);

        // Test hyperbolicCosecant
        double result10 = hyperbolicCosecant(1);
        assert(std::abs(result10 - (1.0 / std::sinh(1))) < 1e-9);

        std::cout << "All unit tests passed." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Unit test failed: " << e.what() << std::endl;
    }
}

int main() {
    try {
        // Example usage of hyperbolic functions
        std::vector<double> values = {0, 1, -1, 2, -2};

        std::cout << "Hyperbolic Sine:" << std::endl;
        for (const auto& value : values) {
            std::cout << "sinh(" << value << ") = " << hyperbolicSine(value) << std::endl;
        }

        std::cout << "Hyperbolic Cosine:" << std::endl;
        for (const auto& value : values) {
            std::cout << "cosh(" << value << ") = " << hyperbolicCosine(value) << std::endl;
        }

        std::cout << "Hyperbolic Tangent:" << std::endl;
        for (const auto& value : values) {
            std::cout << "tanh(" << value << ") = " << hyperbolicTangent(value) << std::endl;
        }

        std::cout << "Hyperbolic Cotangent:" << std::endl;
        for (const auto& value : values) {
            if (value != 0) {
                std::cout << "coth(" << value << ") = " << hyperbolicCotangent(value) << std::endl;
            }
        }

        std::cout << "Hyperbolic Secant:" << std::endl;
        for (const auto& value : values) {
            std::cout << "sech(" << value << ") = " << hyperbolicSecant(value) << std::endl;
        }

        std::cout << "Hyperbolic Cosecant:" << std::endl;
        for (const auto& value : values) {
            if (value != 0) {
                std::cout << "csch(" << value << ") = " << hyperbolicCosecant(value) << std::endl;
            }
        }

        // Example usage of hyperbolic functions for complex numbers
        std::vector<std::complex<double>> complexValues = {std::complex<double>(0, 0), std::complex<double>(1, 1), std::complex<double>(-1, -1)};

        std::cout << "Hyperbolic Sine (Complex):" << std::endl;
        for (const auto& value : complexValues) {
            std::complex<double> result = hyperbolicSineComplex(value);
            displayComplexNumber(result);
        }

        std::cout << "Hyperbolic Cosine (Complex):" << std::endl;
        for (const auto& value : complexValues) {
            std::complex<double> result = hyperbolicCosineComplex(value);
            displayComplexNumber(result);
        }

        std::cout << "Hyperbolic Tangent (Complex):" << std::endl;
        for (const auto& value : complexValues) {
            std::complex<double> result = hyperbolicTangentComplex(value);
            displayComplexNumber(result);
        }

        std::cout << "Hyperbolic Cotangent (Complex):" << std::endl;
        for (const auto& value : complexValues) {
            if (value != std::complex<double>(0, 0)) {
                std::complex<double> result = hyperbolicCotangentComplex(value);
                displayComplexNumber(result);
            }
        }

        std::cout << "Hyperbolic Secant (Complex):" << std::endl;
        for (const auto& value : complexValues) {
            std::complex<double> result = hyperbolicSecantComplex(value);
            displayComplexNumber(result);
        }

        std::cout << "Hyperbolic Cosecant (Complex):" << std::endl;
        for (const auto& value : complexValues) {
            if (value != std::complex<double>(0, 0)) {
                std::complex<double> result = hyperbolicCosecantComplex(value);
                displayComplexNumber(result);
            }
        }

        runUnitTests();
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Unexpected error: " << e.what() << std::endl;
    }

    return 0;
}

#endif // HYPERBOLIC_GEOMETRY_H
