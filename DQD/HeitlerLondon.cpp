#include <functional>

#include "HeitlerLondon.h"

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
                    [setting](double x, double y) {
                        return fockDarwin(x, y, setting, Orientation::Left);
                    });

    SPState right =
            hilbertSpace.createSingleParticleState(
                    [setting](double x, double y) {
                        return fockDarwin(x, y, setting, Orientation::Right);
                    });

    State state_FD_sym = (left ^ right) * sqrt(0.5) +
                         (right ^ left) * sqrt(0.5);

    State state_FD_antisym = (left ^ right) * sqrt(0.5) +
                             (right ^ left) * sqrt(0.5) * (-1.0);

    // Build Hamiltonian

    Operator kineticLeft =
            hilbertSpace.createOperator(
                    [setting](ScalarField field) {
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
                    [setting](ScalarField field) {
                        return kineticEnergy(field, setting);
                    });

    Operator potentialLeft =
            hilbertSpace.createOperator(
                    [setting](ScalarField field) {
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
                    [setting](ScalarField field) {
                        return potentialEnergy(field, setting);
                    });

    Operator coulomb =
            hilbertSpace.createOperator(
                    [setting](double x1, double y1, double x2, double y2) {
                        return coulombEnergy(x1, y1, x2, y2, setting);
                    });

    Operator hamiltonian = kineticLeft + kineticRight + potentialLeft + potentialRight + coulomb;

    // Evaluate energy.

    double energy_sym = hilbertSpace.expectationValue(state_FD_sym, hamiltonian);
    double energy_antisym = hilbertSpace.expectationValue(state_FD_antisym, hamiltonian);

    return (energy_antisym - energy_sym);
}

Complex
fockDarwin(double x, double y, const Setting &setting, Orientation direction) {
    //TODO:
    return Complex();
}

ScalarField
kineticEnergy(ScalarField field, const Setting &setting) {

    return ScalarField(0, 0, 0.0, std::vector<Complex>());
}

ScalarField
potentialEnergy(ScalarField field, const Setting &setting) {
    //TODO:
    return ScalarField(0, 0, 0.0, std::vector<Complex>());
}

Complex
coulombEnergy(double, double, double, double, const Setting &setting) {
    //TODO:
    return Complex();
}
