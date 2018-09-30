#include <cmath>

#include <arrayfire.h>

#include "Fourier.h"


void
fourier::test() {
    af::setDevice(0);
    af::info();
}

std::vector<Complex>
fourier::fft2d(const std::vector<Complex> &field, int width, int height) {
    af::array input(width, height, convertCtoAFC(field).data());
    AFComplex *output = af::fft2(input, width, height).host<AFComplex>();

    return convertAFCtoC(std::vector<AFComplex>(output, output + width * height));
}

std::vector<Complex>
fourier::ifft2d(const std::vector<Complex> &field, int width, int height) {
    af::array input(width, height, convertCtoAFC(field).data());
    AFComplex *output = af::ifft2(input, width, height).host<AFComplex>();

    return convertAFCtoC(std::vector<AFComplex>(output, output + width * height));
}

#include <iostream>

std::vector<Complex>
fourier::convolution(   const std::vector<Complex> & img,    int img_width,     int img_height,
                        const std::vector<Complex> & filter, int filter_width,  int filter_height) {
    af::array input_img(    img_width,      img_height,     convertCtoAFC(img).data());
    af::array input_filter( filter_width,   filter_height,  convertCtoAFC(filter).data());

    int width   = img_width -   filter_width;
    int height  = img_height -  filter_height;

    af::array output = af::fftConvolve2(input_img, input_filter).
            rows(0, width - 1).
            cols(0, height - 1);

    AFComplex *output_data = output.host<AFComplex>();

    return convertAFCtoC(std::vector<AFComplex>(output_data, output_data + width * height));
}

std::vector<AFComplex>
fourier::convertCtoAFC(const std::vector<Complex> &input) {
    std::vector<AFComplex> output(input.size());

    for (int i = 0; i < input.size(); ++i)
        output[i] = AFComplex(input[i].real(), input[i].imag());

    return output;
}

std::vector<Complex>
fourier::convertAFCtoC(const std::vector<AFComplex> &input) {
    std::vector<Complex> output(input.size());

    for (int i = 0; i < input.size(); ++i)
        output[i] = Complex(input[i].real, input[i].imag);

    return output;
}
