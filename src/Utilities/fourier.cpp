#include <cmath>

#include <arrayfire.h>

#include "fourier.hpp"
#include "math_utilities.hpp"


void
fourier::info() {
    af::setDevice(0);
    af::info();
}

ScalarField
fourier::fft2d(const ScalarField &field) {
    int width = field.getWidth();
    int height = field.getHeight();

    af::array input(width, height, convertCtoAFC(field.getDatas()).data());
    AFComplex *output_AFC = af::fft2(input, width, height).host<AFComplex>();

    std::vector<Complex> output_AC = convertAFCtoC(std::vector<AFComplex>(output_AFC, output_AFC + width * height));

    return ScalarField(width, height, field.getGridSize(), output_AC);
}

ScalarField
fourier::ifft2d(const ScalarField &field) {
    int width = field.getWidth();
    int height = field.getHeight();

    af::array input(width, height, convertCtoAFC(field.getDatas()).data());
    AFComplex *output_AFC = af::ifft2(input, width, height).host<AFComplex>();

    std::vector<Complex> output_AC = convertAFCtoC(std::vector<AFComplex>(output_AFC, output_AFC + width * height));

    return ScalarField(width, height, field.getGridSize(), output_AC);
}

ScalarField
fourier::convolution(   const ScalarField & img,
                        const ScalarField & filter) {
    af::array input_img(    img.getWidth(),      img.getHeight(),     convertCtoAFC(img.getDatas()).data());
    af::array input_filter( filter.getWidth(),   filter.getHeight(),  convertCtoAFC(filter.getDatas()).data());

    int img_width   = input_img.dims(0);
    int img_height  = input_img.dims(1);

    int filter_width   = input_filter.dims(0);
    int filter_height  = input_filter.dims(1);

    af::array output = af::fftConvolve2(input_img, input_filter).
            rows(filter_width / 2 - 1 , img_width - filter_width / 2 - 2).
            cols(filter_height / 2 - 1 , img_height - filter_height / 2 - 2);

    int width   = output.dims(0);
    int height  = output.dims(1);

    AFComplex *output_AFC = output.host<AFComplex>();

    std::vector<Complex> output_AC = convertAFCtoC(std::vector<AFComplex>(output_AFC, output_AFC + width * height));

    return reverse(ScalarField(width, height, img.getGridSize(), output_AC));
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
