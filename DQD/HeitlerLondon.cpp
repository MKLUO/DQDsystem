#include <functional>

#include "HeitlerLondon.h"

// TODO: unit tests

// Abbreviations

typedef HilbertSpace::SingleParticleState SPState;
typedef HilbertSpace::State State;

typedef HilbertSpace::Operator Operator;

double
calculateJWithSetting_HL(const Setting &setting) {

    // TODO: build system according to setting

    // TODO: check if setting matches constraint

    HilbertSpace hilbertSpace = HilbertSpace(100, 50, 0.1);

    // Build H-L singlet/triplet state

    SPState left =
            hilbertSpace.createSingleParticleState(
                    fockDarwin(setting, Orientation::Left));

    SPState right =
            hilbertSpace.createSingleParticleState(
                    fockDarwin(setting, Orientation::Right));

    State state_FD_sym = (left ^ right) * sqrt(0.5) +
                         (right ^ left) * sqrt(0.5);

    State state_FD_antisym = (left ^ right) * sqrt(0.5) -
                             (right ^ left) * sqrt(0.5);

    // Build Hamiltonian

    Operator kineticLeft =
            hilbertSpace.createOperator(
                    kineticEnergy(setting),
                    identity());

    Operator kineticRight =
            hilbertSpace.createOperator(
                    identity(),
                    kineticEnergy(setting));


    Operator potentialLeft =
            hilbertSpace.createOperator(
                    potentialEnergy(setting),
                    identity());

    Operator potentialRight =
            hilbertSpace.createOperator(
                    identity(),
                    potentialEnergy(setting));

    Operator coulomb =
            hilbertSpace.createOperator(
                    coulombEnergy(setting));

    Operator hamiltonian =  kineticLeft +
                            kineticRight +
                            potentialLeft +
                            potentialRight +
                            coulomb;

    // Evaluate energy.

    double energy_sym     = hilbertSpace.expectationValue(
                                                state_FD_sym,
                                                hamiltonian);
    double energy_antisym = hilbertSpace.expectationValue(
                                                state_FD_antisym,
                                                hamiltonian);

    return (energy_antisym - energy_sym);
}

SingleParticleScalarFunction
fockDarwin(const Setting &setting, Orientation direction) {
    //TODO:
    return [setting, direction](double x, double y) {
        return Complex();
    };
}

SingleParticleFunction
identity() {
    //TODO: Constants
    return [](ScalarField field) {
        return field;
    };
}


SingleParticleFunction
kineticEnergy(const Setting &setting) {
    //TODO: Constants
    return [setting](ScalarField field) {
        return
            laplacian(field) +
            angularMomentum(field) +
            field * sho_field;
    };

}

SingleParticleFunction
potentialEnergy(const Setting &setting) {
    //TODO: Constants
    return [setting](ScalarField field) {
        return
                field * x_field;
    };
}

DoubleParticleScalarFunction
coulombEnergy(const Setting &setting) {
    //TODO: Constants
    return [setting](double x1, double y1, double x2, double y2) {
        return rInv_field(x1, y1, x2, y2);
    };
}
