#include <iostream>

#include "HilbertSpace.h"
#include "HeitlerLondon.h"
#include "Plot.h"
#include "Fourier.h"

int main() {

    fourier::test();

    int width = 100;
    int height = 40;

    //Create Hilbert Space
    HilbertSpace hilbertSpace = HilbertSpace(width, height, 0.1);

    // Build state
    SPState spState1 =
            hilbertSpace.createSingleParticleState(
                    gaussian(.5));

    State dpState = spState1 ^spState1;

    // Inner product test
    Complex c = dpState * dpState;

    std::cout << c << std::endl;

    // FFT/iFFT test
//    ScalarField field(width, height, 0.1, gaussian(.5));
//
//    std::vector<Complex> field_data = field.getDatas();
//    std::vector<Complex> field_data_FT =
//            fourier::fft2d(field_data, field.getWidth(), field.getHeight());
//    std::vector<Complex> field_data_FT_IFT =
//            fourier::ifft2d(field_data_FT, field.getWidth(), field.getHeight());
//
//    ScalarField field_FT = hilbertSpace.createScalarField(field_data_FT);
//    ScalarField field_FT_IFT = hilbertSpace.createScalarField(field_data_FT_IFT);
//
//    plotter::outputToFile(field, ".\\FIELDCAR");
//    plotter::outputToFile(field_FT, ".\\FIELDCAR2");
//    plotter::outputToFile(field_FT_IFT, ".\\FIELDCAR3");
//
//    // AF array test
//    /*Operator coulomb =
//            hilbertSpace.createOperator(
//                    coulombEnergy(Setting::defaultSetting()));
//
//    double energy_test = hilbertSpace.expectationValue(
//            dpState,
//            coulomb);*/
//
    // Convolution test
    ScalarField img = ScalarField(40, 40, 0.1, scalar(1.));
    ScalarField filter = ScalarField(20, 20, 0.1, scalar(1.));

    ScalarField result = fourier::convolution(img, filter);

    plotter::outputToFile(result, "./CONVCAR");

    return 0;
}
