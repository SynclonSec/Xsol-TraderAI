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

#endif // HYPERBOLIC_GEOMETRY_H

