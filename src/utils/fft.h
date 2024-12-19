
#pragma once

#include <complex>
#include <vector>
#include <cmath>

namespace xsol {
namespace math {

class FFTAnalyzer {
public:
    using Complex = std::complex<double>;
    
    FFTAnalyzer(size_t windowSize = 1024);
    
    // Perform FFT on time series data
    std::vector<Complex> computeFFT(const std::vector<double>& timeSeries);
    
    // Compute power spectrum
    std::vector<double> computePowerSpectrum(const std::vector<Complex>& fftResult);
    
    // Estimate volatility using FFT
    double estimateVolatility(const std::vector<double>& prices, 
                            size_t windowSize = 0);

private:
    size_t m_windowSize;
    std::vector<Complex> m_twiddle; // Pre-computed twiddle factors
    
    void computeTwiddleFactors();
    void bitReverse(std::vector<Complex>& data);
    std::vector<Complex> performFFT(std::vector<Complex>& data);
};

} // namespace math
} // namespace xsol