//#include <arrayfire.h>

#include "Fourier.h"

void
fourier::test() {
    //af::setDevice(0);
    //af::info();
}

std::vector<Complex>
fourier::fft2d(const std::vector<Complex> &, int, int) {
    return std::vector<Complex>();
}

std::vector<Complex>
fourier::ifft2d(const std::vector<Complex> &, int, int) {
    return std::vector<Complex>();
}
