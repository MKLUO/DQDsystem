#include <functional>

#include "HeitlerLondon.h"
#include "MathUtilities.h"

// TODO: unit tests

// Abbreviations

typedef HilbertSpace::SingleParticleState SPState;
typedef HilbertSpace::State State;

typedef HilbertSpace::Operator Operator;

Setting
Setting::defaultSetting() {
    Setting setting;

    setting.width        = 100;
    setting.height       = 50;
    setting.gridSize     = 0.1;
    setting.effectiveE   = 0.;
    setting.effectiveA   = 2;
    setting.B            = 1.;
    setting.alpha        = 1.;
    setting.kappa        = 1.;

    return setting;
}

double
Setting::omegaL() const {
    return Physics::e * B / (2. * Physics::m);
}

double
Setting::omega0() const {
    return alpha * omegaL();
}


double
Setting::omega() const {
    return hypot(omegaL(), omega0());
}

double
Setting::coulombConstant() const {
    return pow(Physics::e, 2.) / (4. * M_PI * Physics::epsilon * Physics::hBar *
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

    Operator hamiltonian =  kineticLeft +
                            kineticRight +
                            potentialLeft +
                            potentialRight +
                            coulomb;

    // TODO: Zeeman energy

    // Evaluate energy.

    double energy_sym     = hilbertSpace.expectationValue(
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
    //TODO:
    return [setting, direction](double x, double y) {
        double a = setting.effectiveA;
        switch (direction) {
            case Orientation::Left:
                return setting.FDConstant() *
                   exp(-0.5i * y * a) *
                   exp(-0.25 * sho_field((x + a), y) * (setting.omega() / setting.omegaL()));
            case Orientation::Right:
                return setting.FDConstant() *
                   exp(+0.5i * y * a) *
                   exp(-0.25 * sho_field((x - a), y) * (setting.omega() / setting.omegaL()));
            default:
                return 0.+0.i;
        }
    };
}

SingleParticleFunction
kineticEnergy(const Setting &setting) {
    return [setting](ScalarField field) {
        return
                laplacian       * field  * -1.0 +
                angularMomentum * field  * 1.i +
                sho_field       * field  * 0.25;
    };
}

SingleParticleFunction
potentialEnergy(const Setting &setting) {
    return [setting](ScalarField field) {
        return
                x_field                     * field * setting.effectiveE +
                quartic(setting.effectiveA) * field * 0.25 * pow(setting.alpha, 2.0);
    };
}

DoubleParticleScalarFunction
coulombEnergy(const Setting &setting) {
    //TODO: Constants
    return [setting](double x1, double y1, double x2, double y2) {
        return
                rInv_field(x1, y1, x2, y2) * setting.coulombConstant();
    };
}
