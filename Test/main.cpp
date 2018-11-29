#include "HL_HM.h"
#include "HilbertSpace.h"
#include "Fourier.h"
#include "Plot.h"

#include <iostream>
#include <corecrt_math_defines.h>

void test1();
void test2();
void test_FD();
void test_CONV();

int main() {
    
    // test1();
    // test2();
    test_FD();
    // test_CONV();

    return 0;
}

void test1() {
    int cases = 20;
    double max_B = 8.;

    std::vector<double> Bs;

    Setting setting = Setting::defaultSetting();
    for (int i = 0; i < cases; ++i)
        Bs.push_back(max_B / double(cases) * double(i));

    for (double B : Bs) {
        setting.B = B;
        std::cout << B << " : " << calculateJWithSetting_HL(setting) << std::endl << std::endl;
    }
}

void test2() {
    
    Setting setting = Setting::defaultSetting();

    HilbertSpace hilbertSpace = HilbertSpace(setting.width, setting.height, setting.gridSize);

    SPState scalarState =
        hilbertSpace.createSingleParticleState(scalar(1.));

    SPState gauss1 =
        hilbertSpace.createSingleParticleState(gaussian(10E-9));

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
        2. * M_PI * 1. / 200E-9,
        2. * M_PI * 0. / 2.));

    ScalarField wave1_lap = laplacian(wave1);

    plotter::outputToFile(wave1, "./FIELDCAR_wave1");
    plotter::outputToFile(wave1_lap, "./FIELDCAR_wave1_lap");
}

void test_FD() {
    Setting setting = Setting::defaultSetting();
    setting.height = 100;
    setting.gridSize = 1.E-9;

    setting.a = 15.E-9;

    HilbertSpace hilbertSpace = HilbertSpace(setting.width, setting.height, setting.gridSize);

    SPState left =
            hilbertSpace.createSingleParticleState(
                    fockDarwin(setting, Orientation::Left));

    SPState right =
            hilbertSpace.createSingleParticleState(
                    fockDarwin(setting, Orientation::Right));

    ScalarField left_sf = left.getField();
    ScalarField right_sf = right.getField();

    plotter::outputToFile(left_sf, "./temp/FIELDCAR_FDL");
    plotter::outputToFile(right_sf, "./temp/FIELDCAR_FDR");
    plotter::outputToFile(left_sf ^ right_sf, "./temp/FIELDCAR_FDLR");         

    Complex test1 = (left ^ right) * (left ^ right);

    Operator i_ts =
            hilbertSpace.createOperator(
                    identity_twoSite(setting));   

    
    Complex result = (i_ts.getOperator()[0])->operatorValue(left ^ right, left ^ right);   

    return;                        
}

void test_CONV() {
    ScalarField c1 = ScalarField(130, 150, 1.0, scalar(1.0));
    ScalarField c2 = ScalarField(100, 100, 1.0, scalar(1.0));

    ScalarField c3 = fourier::convolution(c1, c2);

    plotter::outputToFile(c3, "./temp/C3");
}