#include <math.h>
#include <algorithm>
#include <tuple>

#include "hmp3d.hpp"
#include "math_utilities.hpp"

// TODO: unit tests

#include <iostream>

HMP3D::Setting::Setting(
    const SystemScale scale_, 
    const int basis_, 
    const std::string path_):
    scale(scale_),
    excitedBasis(basis_),
    infoPath(path_) {}

double
HMP3D::Setting::coulombConstant() const {
    return pow(PhysicsContant::e, 2) / (PhysicsContant::kappa * 4. * M_PI * PhysicsContant::epsilon);
}

HMP3D::EigenState::EigenState(
    const ScalarField field_,
    const double energy_) :
    field(field_),
    energy(energy_) {}

std::pair<
    std::vector<HMP3D::EigenState>,
    std::vector<HMP3D::EigenState>>
HMP3D::parseNN_States(std::string infoPath) {
    // TODO: Classify / Sort left & right states
    std::vector<EigenState> states_left, states_right;

    return std::make_pair(states_left, states_right);
}

ScalarField
HMP3D::parseNN_WaveFunction(std::istream is) {
    // TODO: 3D ScalarField
    return ScalarField();
}

Matrix
HMP3D::hamiltonianMatrix(const Setting & setting) {

    HilbertSpace hilbertSpace = HilbertSpace(setting.scale);

    Operator coulomb =
        hilbertSpace.createOperator(
                    coulombEnergy3D(setting));

    std::vector<HMP3D::EigenState> states_lefty, states_right;

    std::tie(states_lefty, states_right) = HMP3D::parseNN_States(setting.infoPath);

    int numBasis = std::min(states_lefty.size(), states_right.size());
    if (setting.excitedBasis < numBasis - 1)
        numBasis = setting.excitedBasis + 1;

    
    // TODO: Add detuning and Zeeman.
    // TODO: Check energies of (1, 1) vs. (2, 0)

    std::vector<State> states;
    std::vector<double> energies;
    
    // (1, 1)
    for (int i = 0; i < numBasis; ++i)
    for (Spin::Type spin_lefty : {Spin::Type::Up, Spin::Type::Down})
    for (int j = 0; j < numBasis; ++j)
    for (Spin::Type spin_right : {Spin::Type::Up, Spin::Type::Down}) {

        states.push_back(State(
            SPState(states_lefty[i].field, Spin(spin_lefty), "L" + std::to_string(i)),
            SPState(states_right[j].field, Spin(spin_right), "R" + std::to_string(j))).antisym());

        energies.push_back(
            states_lefty[i].energy + 
            states_right[j].energy);
    }

    // (2, 0) & (0, 2)
    states.push_back(State(
        SPState(states_lefty[0].field, Spin(Spin::Type::Up), "L0"),
        SPState(states_lefty[0].field, Spin(Spin::Type::Down), "L0")).antisym());
    energies.push_back(states_lefty[0].energy * 2.0);

    states.push_back(State(
        SPState(states_right[0].field, Spin(Spin::Type::Up), "R0"),
        SPState(states_right[0].field, Spin(Spin::Type::Down), "R0")).antisym());
    energies.push_back(states_right[0].energy * 2.0);

    int numStates = states.size();

    Matrix matrix;
    for (int i = 0; i < numStates; ++i) {
        std::vector<ComplexHighRes> vec;
        for (int j = 0; j < numStates; ++j) {
            ComplexHighRes energy = hilbertSpace.operatorValue(states[i], coulomb, states[j]);
            if (i == j) energy += Complex(energies[i]);
        }
        matrix.push_back(vec);
    }

    return matrix;
}