#include "HeitlerLondon.h"
#include "HilbertSpace.h"

#include "Plot.h"

#include <iostream>

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

    plotter::outputToFile(gauss_lap, "./FIELDCAR_FIELDCAR_gausslap");
    plotter::outputToFile(gauss1.getField(), "./FIELDCAR_gauss");

    Complex a = hilbertSpace.operatorValue(state2, identityOp + identityOp, state1);
    Complex b = hilbertSpace.expectationValue(state2, identityOp + identityOp);
    Complex c = hilbertSpace.expectationValue(state3, identityOp);
    Complex d = hilbertSpace.expectationValue(state3, laplaceOp1);

    return;
}