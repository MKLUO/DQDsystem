#include "HeitlerLondon.h"
#include "HilbertSpace.h"

#include "Plot.h"

#include <iostream>
#include <corecrt_math_defines.h>

void test1();
void test2();

int main() {
    
    //test1();
    test2();

    return 0;
}

void test1() {
    int cases = 5;
    double max_B = 8.;

    std::vector<double> Bs;

    Setting setting = Setting::defaultSetting();
    for (int i = 0; i < cases; ++i)
        Bs.push_back(max_B / double(cases) * double(i));

    for (double B : Bs) {
        setting.B = B;
        std::cout << calculateJWithSetting_HL(setting) << std::endl;
    }
}

void test2() {
    
    Setting setting = Setting::defaultSetting();
    setting.width = 100;
    setting.height = 100;
    setting.gridSize = 0.02;

    HilbertSpace hilbertSpace = HilbertSpace(setting.width, setting.height, setting.gridSize);

    SPState scalarState =
        hilbertSpace.createSingleParticleState(scalar(1.));

    SPState gauss1 =
        hilbertSpace.createSingleParticleState(gaussian(0.2));

    State state1 = (scalarState ^ scalarState).normalize();
    State state2 = (scalarState ^ scalarState).normalize() * 2;

    State state3 = (gauss1 ^ scalarState);


    Operator identityOp =
        hilbertSpace.createOperator(
            identity,
            identity);

    Operator laplaceOp1 =
        hilbertSpace.createOperator(
            laplacian,
            identity);

    ScalarField gauss_lap = laplacian(gauss1.getField());

    Complex a = hilbertSpace.operatorValue(state2, identityOp + identityOp, state1);
    Complex b = hilbertSpace.expectationValue(state2, identityOp + identityOp);
    Complex c = hilbertSpace.expectationValue(state3, identityOp);
    Complex d = hilbertSpace.expectationValue(state3, laplaceOp1);

    plotter::outputToFile(gauss_lap, "./FIELDCAR_gausslap");
    plotter::outputToFile(gauss1.getField(), "./FIELDCAR_gauss");

    ScalarField wave1 = hilbertSpace.createScalarField(planeWave(   
        2. * M_PI * 80. / 2.,
        2. * M_PI * 10. / 2.));

    ScalarField wave1_lap = laplacian(wave1);

    plotter::outputToFile(wave1, "./FIELDCAR_wave1");
    plotter::outputToFile(wave1_lap, "./FIELDCAR_wave1_lap");

    return;
}