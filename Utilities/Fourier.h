#pragma once

#include "MathUtilities.h"

namespace fourier {
    void test();

    std::vector<Complex>
    fft2d(const std::vector<Complex> &, int, int);

    std::vector<Complex>
    ifft2d(const std::vector<Complex> &, int, int);
}