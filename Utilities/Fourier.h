#pragma once

#include <arrayfire.h>

#include "MathUtilities.h"

using AFComplex = af::cdouble;

namespace fourier {
    void info();

    ScalarField
    fft2d(const ScalarField &);

    ScalarField
    ifft2d(const ScalarField &);

    ScalarField
    convolution(const ScalarField &, const ScalarField &);

    std::vector<AFComplex>
    convertCtoAFC(const std::vector<Complex> &);

    std::vector<Complex>
    convertAFCtoC(const std::vector<AFComplex> &);
}