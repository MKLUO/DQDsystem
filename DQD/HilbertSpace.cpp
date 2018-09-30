#include <utility>

#include "HilbertSpace.h"

//////////////////////////////
//       HilbertSpace       //
//////////////////////////////

HilbertSpace::HilbertSpace(int width_, int height_, double gridSize_) {
    width = width_;
    height = height_;
    gridSize = gridSize_;
}

SPState
HilbertSpace::createSingleParticleState(const SPSFunction &function) const {
    return SPState(ScalarField(width, height, gridSize, function));
}

Operator
HilbertSpace::createOperator(const SingleParticleFunction &function1, const SingleParticleFunction &function2) const {
    return Operator({new SPOperator(function1, function2)});
}

Operator
HilbertSpace::createOperator(const DoubleParticleScalarFunction &function) const {
    return Operator({new DPSOperator(function)});
}

Complex
HilbertSpace::operatorValue(const State &stateLeft, const Operator &ops, const State &stateRight) const {

    Complex result;

    for (SingleOperator *op : ops.getOperator()) {
        result += op->operatorValue(stateLeft, stateRight);
    }

    return result;
}

double
HilbertSpace::expectationValue(const State &state, const Operator &ops) const {
    return HilbertSpace::operatorValue(state, ops, state).real();
}

//////////////////////////////
//   SingleParticleState    //
//////////////////////////////
SPState::SingleParticleState(const ScalarField &field_) :
        field(field_) {}

ScalarField
SPState::getField() const {
    return field;
}

SPState
SPState::operator+(const SPState &state) const {
    return SPState(this->getField() + state.getField());
}

SPState
SPState::operator*(Complex c) const {
    return SPState(field * c);
}

Complex
SPState::operator*(const SPState &state) const {
    return (this->getField() * state.getField());
}

State
SPState::operator^(const SPState &state) const {
    return State(*this, state);
}

//////////////////////////////
// SingleParticleStatePair  //
//////////////////////////////

SPStatePair::SingleParticleStatePair(const SPState &state1, const SPState &state2) :
        first(state1),
        second(state2) {}

SPState
SPStatePair::getFirstField() const {
    return first;
}

SPState
SPStatePair::getSecondField() const {
    return second;
}

SPStatePair
SPStatePair::operator*(Complex c) const {
    return SPStatePair(first * c, second);
}


//////////////////////////////
//          State           //
//////////////////////////////

State::State(const SPState &state1, const SPState &state2) :
        states({SingleParticleStatePair(state1, state2)}) {}

State::State(const std::vector<SPStatePair> &states_) {
    states = states_;
}

std::vector<SPStatePair>
State::getState() const {
    return states;
}

State
State::operator+(const State &state) const {
    std::vector<SPStatePair> v1 = this->getState();
    std::vector<SPStatePair> v2 = state.getState();

    v1.insert(v1.end(), v2.begin(), v2.end());

    return State(v1);
}

State
State::operator-(const State &state) const {
    std::vector<SPStatePair> v1 = this->getState();
    std::vector<SPStatePair> v2 = state.getState();

    for (SPStatePair &pair : v2)
        pair = pair * Complex(-1., 0.);

    v1.insert(v1.end(), v2.begin(), v2.end());

    return State(v1);
}

State
State::operator*(Complex c) const {
    std::vector<SPStatePair> states = this->getState();
    for (SPStatePair &pair : states)
        pair = pair * c;
    return State(states);
}

Complex
State::operator*(const State &state) const {
    std::vector<SPStatePair> states1 = this->getState();
    std::vector<SPStatePair> states2 = state.getState();

    Complex result;

    for (const SPStatePair &pair1 : states1) {
        for (const SPStatePair &pair2 : states2) {
            SPState field1Left = pair1.getFirstField();
            SPState field1Right = pair1.getSecondField();
            SPState field2Left = pair2.getFirstField();
            SPState field2Right = pair2.getSecondField();

            result += (field1Left * field2Left) * (field1Right * field2Right);
        }
    }

    return result;
}

//////////////////////////////
//  SingleParticleOperator  //
//////////////////////////////

SPOperator::SingleParticleOperator(const SingleParticleFunction &left_,
                                   const SingleParticleFunction &right_) {
    left = left_;
    right = right_;
}

State
SPOperator::operator*(const State &state) const {
    std::vector<SPStatePair> states = state.getState();

    for (SPStatePair &pair : states)
        pair = *this * pair;

    return State(states);
}

SPStatePair
SPOperator::operator*(const SPStatePair &pair) const {
    return SPStatePair(SPState(pair.getFirstField().getField() * left),
                       SPState(pair.getSecondField().getField() * right));
}

Complex
SPOperator::operatorValue(const State &left, const State &right) const {
    return left * (*this * right);
}

//////////////////////////////////
// DoubleParticleScalarOperator //
//////////////////////////////////

DPSOperator::DoubleParticleScalarOperator(const DoubleParticleScalarFunction &func_) {
    func = func_;
}

/*
State
DPSOperator::operator*(const State &) const {
    return State(std::vector());
}
*/

Complex
DPSOperator::operatorValue(const State &left, const State &right) const {
    std::vector<SPStatePair> states1 = left.getState();
    std::vector<SPStatePair> states2 = right.getState();

    Complex result;

    for (const SPStatePair &pair1 : states1) {
        for (const SPStatePair &pair2 : states2) {
            ScalarField field1Left = pair1.getFirstField().getField();
            ScalarField field1Right = pair1.getSecondField().getField();
            ScalarField field2Left = pair2.getFirstField().getField();
            ScalarField field2Right = pair2.getSecondField().getField();

            result += twoSiteIntegral(field1Left, field2Left, func, field1Right, field2Right);
        }
    }

    return result;
}

//////////////////////////////
//         Operator         //
//////////////////////////////

Operator::Operator(const std::vector<SingleOperator *> &operators_) {
    operators = operators_;
}

Operator::~Operator() {
    for (auto op : operators)
        delete (op);
}

Operator
Operator::operator+(const Operator &ops) const {
    std::vector<SingleOperator *> v1 = this->getOperator();
    std::vector<SingleOperator *> v2 = ops.getOperator();

    v1.insert(v1.end(), v2.begin(), v2.end());

    return Operator(v1);
}

std::vector<SingleOperator *>
Operator::getOperator() const {
    return operators;
}