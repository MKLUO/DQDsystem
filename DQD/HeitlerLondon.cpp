#include <functional>

#include "HeitlerLondon.h"

// Abbreviations

typedef HilbertSpace::SingleParticleState SPState;
typedef HilbertSpace::SingleParticleStatePair SPStatePair;
typedef HilbertSpace::State State;

typedef HilbertSpace::Operator Operator;

double calculateJWithSetting_HL(const Setting &setting) {

    HilbertSpace hilbertSpace = HilbertSpace(100, 50, 0.1);

    // Build H-L singlet/triplet state

    SPState left =
            hilbertSpace.createSingleParticleState(
                    [](int x, int y) {
                        return fockDarwin(x, y, setting, Orientation::Left);
                    });

    SPState right =
            hilbertSpace.createSingleParticleState(
                    [](int x, int y) {
                        return fockDarwin(x, y, setting, Orientation::Right);
                    });

    State state_FD_sym      =   State({left ^ right}) * sqrt(0.5) +
                                State({right ^ left}) * sqrt(0.5);

    State state_FD_antisym  =   State({left ^ right}) * sqrt(0.5) +
                                State({right ^ left}) * sqrt(0.5) * (-1.0);

    // Build Hamiltonian

    Operator kineticLeft =
            hilbertSpace.createOperator(
                    [](ScalarField field) {
                        return kineticEnergy(field, setting);
                    },
                    [](ScalarField field) {
                        return field;
                    });

    Operator kineticRight =
            hilbertSpace.createOperator(
                    [](ScalarField field) {
                        return field;
                    },
                    [](ScalarField field) {
                        return kineticEnergy(field, setting);
                    });

    Operator potentialLeft =
            hilbertSpace.createOperator(
                    [](ScalarField field) {
                        return potentialEnergy(field, setting);
                    },
                    [](ScalarField field) {
                        return field;
                    });

    Operator potentialRight =
            hilbertSpace.createOperator(
                    [](ScalarField field) {
                        return field;
                    },
                    [](ScalarField field) {
                        return potentialEnergy(field, setting);
                    });

    Operator coulomb =
            hilbertSpace.createOperator(
                    [](int x1, int y1, int x2, int y2) {
                        return coulombEnergy(x1, y1, x2, y2, setting);
                    });

    Operator hamiltonian = kineticLeft + kineticRight + potentialLeft + potentialRight + coulomb;

    // Evaluate energy.

    double energy_sym = hilbertSpace.expectationValue(state_FD_sym, hamiltonian);
    double energy_antisym = hilbertSpace.expectationValue(state_FD_antisym, hamiltonian);

    return (energy_antisym - energy_sym);
}

Complex fockDarwin(int x, int y, const Setting &setting, Orientation direction) {
    return Complex();
}

ScalarField kineticEnergy(ScalarField field, const Setting &setting) {
    return ScalarField();
}