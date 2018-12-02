#include <functional>

#define _USE_MATH_DEFINES

#include <math.h>

#include "hl_hm.hpp"
#include "math_utilities.hpp"

// TODO: unit tests

#include <iostream>

Setting
Setting::defaultSetting() {

    Setting setting;

    setting.scale = HilbertSpace::SystemScale::defaultScale();
    
    setting.a = 20E-9;
    setting.d = 100E-9;
    setting.B = 1.0;

    return setting;
}

double
Setting::omegaL() const {
    return PhysicsContant::e * B / (2. * PhysicsContant::m);
}

double
Setting::omega0() const {
    return PhysicsContant::hBar / (PhysicsContant::m * pow(a, 2));
}

double
Setting::omega() const {
    return hypot(omegaL(), omega0());
}

double
Setting::coulombConstant() const {
    return pow(PhysicsContant::e, 2) / (PhysicsContant::kappa * 4. * M_PI * PhysicsContant::epsilon);
}

double
Setting::FDConstant() const {
    return sqrt(PhysicsContant::m * omega() / (M_PI * PhysicsContant::hBar));
}

double
Setting::magneticLength() const {
    return sqrt(PhysicsContant::hBar / (PhysicsContant::e * B));
}

Matrix
coulombMatrix(const Setting & setting, Ansatz ansatz, Basis basis) {

    HilbertSpace hilbertSpace = HilbertSpace(setting.scale);

    Operator coulomb =
        hilbertSpace.createOperator(
                    coulombEnergy(setting));

    SPState left       = fockDarwinWithSpin(hilbertSpace, setting, Orientation::Left,  Spin::None);
    SPState right      = fockDarwinWithSpin(hilbertSpace, setting, Orientation::Right, Spin::None);
    SPState left_up    = fockDarwinWithSpin(hilbertSpace, setting, Orientation::Left,  Spin::Up);
    SPState left_down  = fockDarwinWithSpin(hilbertSpace, setting, Orientation::Left,  Spin::Down);
    SPState right_up   = fockDarwinWithSpin(hilbertSpace, setting, Orientation::Right, Spin::Up);
    SPState right_down = fockDarwinWithSpin(hilbertSpace, setting, Orientation::Right, Spin::Down);

    std::vector<State> states;

    switch (ansatz)
    {
        case Ansatz::HL:
            states.push_back((left ^ right).sym());
            states.push_back((left ^ right).antisym());
            break;

        case Ansatz::HM:
            states.push_back((left_up   ^ right_down).antisym());
            states.push_back((left_down ^ right_up)  .antisym());
            states.push_back((left_up   ^ left_down) .antisym());
            states.push_back((right_up  ^ right_down).antisym());
            break;
    }                    

    Matrix matrix;
    for (const State& state1 : states) {
        std::vector<ComplexHighRes> vec;
        for (const State& state2 : states)
            vec.push_back(hilbertSpace.operatorValue(state1, coulomb, state2));
        matrix.push_back(vec);
    }

    return matrix;
}

double
calculateJWithSetting_HL(const Setting &setting) {

    // TODO: check if setting matches constraint
    // TODO: check if system is too small that the F-D wave functions are not covered.

    HilbertSpace hilbertSpace = HilbertSpace(setting.scale);

    // Build H-L singlet/triplet state

    SPState left =
            hilbertSpace.createSingleParticleState(
                    fockDarwin(setting, Orientation::Left));

    SPState right =
            hilbertSpace.createSingleParticleState(
                    fockDarwin(setting, Orientation::Right));

    // TODO: Normalization!
    State state_FD_sym =        ((left ^ right) + (right ^ left)).normalize();

    State state_FD_antisym =    ((left ^ right) - (right ^ left)).normalize();

    // Build Hamiltonian

    Operator kineticLeft =
            hilbertSpace.createOperator(
                    kineticEnergy(setting),
                    identity);

    Operator kineticRight =
            hilbertSpace.createOperator(
                    identity,
                    kineticEnergy(setting));


    Operator potentialLeft =
            hilbertSpace.createOperator(
                    potentialEnergy(setting),
                    identity);

    Operator potentialRight =
            hilbertSpace.createOperator(
                    identity,
                    potentialEnergy(setting));

    Operator coulomb =
            hilbertSpace.createOperator(
                    coulombEnergy(setting));

    Operator hamiltonian = kineticLeft +
                           kineticRight +
                           potentialLeft +
                           potentialRight +
                           coulomb;

    // TODO: Zeeman energy

    // Evaluate energy.


    // TEST
    // State state_FD_LR = (left ^ right).normalize();
    // State state_FD_RL = (right ^ left).normalize();
    // Complex test1 = (coulomb.getOperator()[0])->operatorValue(state_FD_LR, state_FD_LR);
    // Complex test2 = (coulomb.getOperator()[0])->operatorValue(state_FD_RL, state_FD_RL);
    // Complex test3 = (coulomb.getOperator()[0])->operatorValue(state_FD_LR, state_FD_RL);
    // Complex test4 = (coulomb.getOperator()[0])->operatorValue(state_FD_RL, state_FD_LR);
    // Complex testss = (coulomb.getOperator()[0])->operatorValue(state_FD_sym, state_FD_sym);
    // Complex testsa = (coulomb.getOperator()[0])->operatorValue(state_FD_antisym, state_FD_antisym);

    // TODO: Float point data distortion!!!

    ComplexHighRes result_sym = 0., result_asym = 0.;

    // for (SingleOperator *op : hamiltonian.getOperator()) {
    //     Complex temp1 = op->operatorValue(state_FD_sym, state_FD_sym).real();
    //     Complex temp2 = op->operatorValue(state_FD_antisym, state_FD_antisym).real();
    //     result_sym += temp1;
    //     result_asym += temp2;

    //     //std::cout << temp1 << " " << temp2 << std::endl;
    // }

    result_sym = hilbertSpace.operatorValue(state_FD_sym, hamiltonian, state_FD_sym);
    result_asym = hilbertSpace.operatorValue(state_FD_antisym, hamiltonian, state_FD_antisym);

    return (result_asym - result_sym).real();
}

// ScalarFunctions

SingleParticleScalarFunction
fockDarwin(const Setting &setting, Orientation direction) {
    return [setting, direction](double x, double y) {
        double d = setting.d;
        double B = setting.B;
        double lb = setting.magneticLength();
        double omega = setting.omega();
        double sign;
        switch (direction) {
            case Orientation::Left:
                sign = -1.0;
                break;
            case Orientation::Right:
                sign = +1.0;
                break;
        }
        if (B == 0)
            return setting.FDConstant() *
                exp(- PhysicsContant::m * omega * sho_field((x - d * sign / 2), y) / 
                    (2.0 * PhysicsContant::hBar));
        else
            return setting.FDConstant() *
                exp(1.i * sign * y * d  / (4.0 * pow(lb, 2))) *
                exp(- PhysicsContant::m * omega * sho_field((x - d * sign / 2), y) / 
                    (2.0 * PhysicsContant::hBar));
    };
}

SingleParticleFunction
kineticEnergy(const Setting &setting) {
    return [setting](ScalarField field) {
        double B = setting.B;
        return  1.0 / (2.0 * PhysicsContant::m) * (
                laplacian * field * (- pow(PhysicsContant::hBar, 2)) +
                angularMomentum * field * (1.i * PhysicsContant::hBar * PhysicsContant::e * B)  +
                sho_field * field * (pow(0.5 * PhysicsContant::e * B, 2)) );
    };
}

SingleParticleFunction
potentialEnergy(const Setting &setting) {
    return [setting](ScalarField field) {
        double d = setting.d;
        double omega0 = setting.omega0();
        return quartic(0.5 * d) * field * (0.5 * PhysicsContant::m * pow(omega0, 2));
    };
}

DoubleParticleScalarFunction
coulombEnergy(const Setting &setting) {
    return [setting](double x1, double y1, double x2, double y2) {
        return
                rInv_field(setting.scale.gridSize)(x1, y1, x2, y2) * setting.coulombConstant();
    };
}

// Util

SPState
fockDarwinWithSpin(const HilbertSpace & hilbertSpace, Setting setting, Orientation orient, Spin spin) {
    return hilbertSpace.createSingleParticleState(
                    fockDarwin(setting, orient), spin);
}


// For test
DoubleParticleScalarFunction
identity_twoSite(const Setting &setting) {
    return [](double x1, double y1, double x2, double y2) {
        return 1.0;
    };
}
