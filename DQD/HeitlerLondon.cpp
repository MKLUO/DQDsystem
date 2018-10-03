#include <functional>

#define _USE_MATH_DEFINES

#include <math.h>

#include "HeitlerLondon.h"
#include "MathUtilities.h"

// TODO: unit tests

Setting
Setting::defaultSetting() {
    Setting setting;

    setting.width = 200;
    setting.height = 100;
    setting.gridSize = 0.01;
    setting.E = 0.0;
    setting.a = 0.3;
    setting.B = 0.000000000001;
    setting.alpha = 10.;
    setting.kappa = 11.7;

    return setting;
}

double
Setting::omegaL() const {
    return Physics::e * B / (2. * Physics::m);
}

double
Setting::omegaL(double B) const {
    return Physics::e * B / (2. * Physics::m);
}

double
Setting::omega0() const {
    return alpha * omegaL(1.0);
}

double
Setting::omega() const {
    return hypot(omegaL(), omega0());
}

double
Setting::coulombConstant() const {
    return pow(B, -0.5) * pow(Physics::e, 2.) / (4. * M_PI * Physics::epsilon * Physics::hBar *
                                  kappa * magneticLength() * omegaL());
}

double
Setting::FDConstant() const {
    return sqrt(Physics::m * omega() / (M_PI * Physics::hBar));
}

double
Setting::magneticLength() const {
    return sqrt(Physics::hBar / (Physics::e * B));
}

double
Setting::magneticLength(double B) const {
    return sqrt(Physics::hBar / (Physics::e * B));
}

double
calculateJWithSetting_HL(const Setting &setting) {

    // TODO: check if setting matches constraint

    HilbertSpace hilbertSpace = HilbertSpace(setting.width, setting.height, setting.gridSize);

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

    double energy_sym = hilbertSpace.expectationValue(
            state_FD_sym,
            hamiltonian);
    double energy_antisym = hilbertSpace.expectationValue(
            state_FD_antisym,
            hamiltonian);

    return (energy_antisym - energy_sym) * Physics::hBar * setting.omegaL();
}

// ScalarFunctions

SingleParticleScalarFunction
fockDarwin(const Setting &setting, Orientation direction) {
    return [setting, direction](double x, double y) {
        double a = setting.a;
        double B = setting.B;
        double sign;
        switch (direction) {
            case Orientation::Left:
                sign = +1.0;
                break;
            case Orientation::Right:
                sign = -1.0;
                break;
        }
        return setting.FDConstant() * setting.magneticLength(1.) *
            exp(Complex(-0.5i) * y * a * sign * B) *
            exp(-0.25 * sho_field((x + a * sign), y) * B * (setting.omega() / setting.omegaL()));
    };
}

SingleParticleFunction
kineticEnergy(const Setting &setting) {
    return [setting](ScalarField field) {
        double B = setting.B;
        return
                laplacian * field * (-1.0 / B) +
                angularMomentum * field * Complex(1.i) +
                sho_field * field * 0.25 * B;
    };
}

SingleParticleFunction
potentialEnergy(const Setting &setting) {
    return [setting](ScalarField field) {
        double a = setting.a;
        double B = setting.B;
        double E = setting.E;
        return
                x_field * field * (E / B) +
                quartic(a) * field * 0.25 * (pow(setting.alpha, 2.0) / B);
    };
}

DoubleParticleScalarFunction
coulombEnergy(const Setting &setting) {
    //TODO: Constants
    return [setting](double x1, double y1, double x2, double y2) {
        return
                rInv_field(setting.gridSize)(x1, y1, x2, y2) * setting.coulombConstant();
    };
}
