#ifndef EULER_H
#define EULER_H

#include <iostream>
#include <vector>
#include <complex>
#include <string>

void log(const std::string& message);
std::complex<double> eulerFormula(double x);
std::vector<std::complex<double>> dft(const std::vector<double>& priceData);
void displayComplexNumber(const std::complex<double>& number);
void runUnitTests();

#endif // EULER_H
