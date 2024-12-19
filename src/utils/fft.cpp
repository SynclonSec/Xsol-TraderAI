

#include "fft.h"
#include <algorithm>
#include <stdexcept>

namespace xsol {
namespace math {

FFTAnalyzer::FFTAnalyzer(size_t windowSize) 
    : m_windowSize(windowSize) {
    // Ensure window size is power of 2
    if (m_windowSize & (m_windowSize - 1)) {
        throw std::invalid_argument("Window size must be power of 2");
    }
    computeTwiddleFactors();
}

void FFTAnalyzer::computeTwiddleFactors() {
    m_twiddle.resize(m_windowSize);
    for (size_t i = 0; i < m_windowSize; ++i) {
        double angle = -2.0 * M_PI * i / m_windowSize;
        m_twiddle[i] = Complex(cos(angle), sin(angle));
    }
}

void FFTAnalyzer::bitReverse(std::vector<Complex>& data) {
    size_t n = data.size();
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        while (j >= bit) {
            j -= bit;
            bit >>= 1;
        }
        j += bit;
        if (i < j) {
            std::swap(data[i], data[j]);
        }
    }
}

std::vector<FFTAnalyzer::Complex> FFTAnalyzer::performFFT(std::vector<Complex>& data) {
    size_t n = data.size();
    bitReverse(data);
    
    // Cooley-Tukey FFT algorithm
    for (size_t len = 2; len <= n; len <<= 1) {
        size_t halfLen = len >> 1;
        for (size_t i = 0; i < n; i += len) {
            for (size_t j = 0; j < halfLen; ++j) {
                Complex temp = data[i + j + halfLen] * 
                             m_twiddle[j * m_windowSize / len];
                data[i + j + halfLen] = data[i + j] - temp;
                data[i + j] = data[i + j] + temp;
            }
        }
    }
    
    return data;
}

std::vector<FFTAnalyzer::Complex> FFTAnalyzer::computeFFT(
    const std::vector<double>& timeSeries) {
    
    if (timeSeries.size() != m_windowSize) {
        throw std::invalid_argument("Input size must match window size");
    }
    
    std::vector<Complex> data(m_windowSize);
    for (size_t i = 0; i < m_windowSize; ++i) {
        data[i] = Complex(timeSeries[i], 0.0);
    }
    
    return performFFT(data);
}

std::vector<double> FFTAnalyzer::computePowerSpectrum(
    const std::vector<Complex>& fftResult) {
    
    std::vector<double> powerSpectrum(fftResult.size() / 2 + 1);
    
    for (size_t i = 0; i < powerSpectrum.size(); ++i) {
        powerSpectrum[i] = std::abs(fftResult[i]) * std::abs(fftResult[i]) 
                          / m_windowSize;
    }
    
    return powerSpectrum;
}

double FFTAnalyzer::estimateVolatility(
    const std::vector<double>& prices, 
    size_t windowSize) {
    
    if (windowSize != 0) {
        m_windowSize = windowSize;
        computeTwiddleFactors();
    }
    
    // Compute returns
    std::vector<double> returns(m_windowSize);
    for (size_t i = 1; i < prices.size() && i < m_windowSize; ++i) {
        returns[i] = log(prices[i] / prices[i-1]);
    }
    
    // Compute FFT
    auto fftResult = computeFFT(returns);
    
    // Compute power spectrum
    auto powerSpectrum = computePowerSpectrum(fftResult);
    
    // Estimate volatility from power spectrum
    double totalPower = 0.0;
    for (size_t i = 1; i < powerSpectrum.size(); ++i) {
        totalPower += powerSpectrum[i] * i; // Weighted by frequency
    }
    
    return sqrt(totalPower * 252.0); // Annualized volatility
}

} // namespace math
} // namespace xsol
