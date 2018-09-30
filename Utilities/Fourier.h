#pragma once

#include <arrayfire.h>

#include "MathUtilities.h"

using AFComplex = af::cdouble;

namespace fourier {
    void test();

    std::vector<Complex>
    fft2d(const std::vector<Complex> &, int, int);

    std::vector<Complex>
    ifft2d(const std::vector<Complex> &, int, int);

    std::vector<Complex>
    convolution(const std::vector<Complex> &, int, int, const std::vector<Complex> &, int, int);

    std::vector<AFComplex>
    convertCtoAFC(const std::vector<Complex> &);

    std::vector<Complex>
    convertAFCtoC(const std::vector<AFComplex> &);
}