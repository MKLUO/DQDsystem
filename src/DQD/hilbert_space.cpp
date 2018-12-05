#include <utility>

#include "hilbert_space.hpp"

//////////////////////////////
//       HilbertSpace       //
//////////////////////////////

HilbertSpace::HilbertSpace(
    SystemScale scale_):
    scale(scale_) {}

SPState
HilbertSpace::createSingleParticleState(const SPSFunction &function, const Spin & spin, const std::string & label) const {
    if (label == "") 
        return SPState(createScalarField(function), spin).normalize();
    else 
        return SPState(createScalarField(function), spin, label).normalize();
}

ScalarField
HilbertSpace::createScalarField() const {
    return ScalarField(scale);
}

ScalarField
HilbertSpace::createScalarField(const std::vector<Complex> &field) const {
    return ScalarField(scale, field);
}

ScalarField
HilbertSpace::createScalarField(const SPSFunction &function) const {
    return ScalarField(scale, function);
}


Operator
HilbertSpace::createOperator(const SingleParticleFunction &function1, const SingleParticleFunction &function2) const {
    return Operator({new SPOperator(function1, function2)});
}

Operator
HilbertSpace::createOperator(const DoubleParticleScalarFunction &function) const {
    return Operator({new DPSOperator(function)});
}

ComplexHighRes
HilbertSpace::operatorValue(const State &stateLeft, const Operator &ops, const State &stateRight) const {

    ComplexHighRes result;

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

int SPState::auto_index = 1;

HilbertSpace::SingleParticleState::SingleParticleState(
    const ScalarField &field_, 
    const Spin & spin_, 
    const std::string & label_):
        field(field_), 
        spin(spin_), 
        label(
            (label_ == "")?
                ("(" + std::to_string(auto_index++) + ")"):
                label_) {}        

ScalarField
HilbertSpace::SingleParticleState::getField() const {
    return field;
}

Spin
HilbertSpace::SingleParticleState::getSpin() const {
    return spin;
}

std::string
HilbertSpace::SingleParticleState::getLabel() const {
    return label;
}

SPState
HilbertSpace::SingleParticleState::normalize() const {
    return SPState(field.normalize(), spin, label);
}

SPState
HilbertSpace::SingleParticleState::operator+(const SPState &state) const {
    if (spin == state.spin)
        return SPState(this->getField() + state.getField(), spin, label + "\'");
    else 
        throw std::exception();
}

SPState
HilbertSpace::SingleParticleState::operator-(const SPState &state) const {
    return operator+(state * (-1.));
}

SPState
HilbertSpace::SingleParticleState::operator*(Complex c) const {
    return SPState(field * c, spin, label);
}

SPState
HilbertSpace::SingleParticleState::operator/(Complex c) const {
    return operator*(1./c);
}

// Spin selection rule is implemented here. TODO: Is it a good practice?
ComplexHighRes
HilbertSpace::SingleParticleState::operator*(const SPState &state) const {
    if ((spin == Spin::None) && (state.spin == Spin::None)) {
        return (this->getField() * state.getField());

    } else if ((spin != Spin::None) && (state.spin != Spin::None)) {
        if (spin == state.spin)
            return (this->getField() * state.getField());            
        else 
            return 0.0;

    } else {
        throw std::exception();
    }
}

State
HilbertSpace::SingleParticleState::operator^(const SPState &state) const {
    return State(*this, state);
}

//////////////////////////////
// SingleParticleStatePair  //
//////////////////////////////

HilbertSpace::SingleParticleStatePair::SingleParticleStatePair(
    const SPState &state1, 
    const SPState &state2, 
    const Complex & coef_,
    const std::string & label) :
        first(state1),
        second(state2),
        coef(coef_),
        label_override(label) {}

SPState
HilbertSpace::SingleParticleStatePair::getFirstField() const {
    return first;
}

SPState
HilbertSpace::SingleParticleStatePair::getSecondField() const {
    return second;
}

Complex
HilbertSpace::SingleParticleStatePair::getCoef() const {
    return coef;
}

std::string
HilbertSpace::SingleParticleStatePair::getLabel() const {
    if (label_override != "") return label_override;
 
    return (
        std::to_string(coef.real()) + 
        "[" + 
        first.getLabel() + spinSign(first.getSpin()) + " " + 
        second.getLabel() + spinSign(second.getSpin()) + "]");
}

SPStatePair
HilbertSpace::SingleParticleStatePair::operator*(Complex c) const {
    return SPStatePair(first, second, coef * c);
}


//////////////////////////////
//          State           //
//////////////////////////////

HilbertSpace::State::State(const SPState &state1, const SPState &state2, const std::string & label) :
        states({SingleParticleStatePair(state1, state2)}), 
        label_override(label) {}

HilbertSpace::State::State(const std::vector<SPStatePair> &states_, const std::string & label) :
        states(states_),
        label_override(label) {}

std::vector<SPStatePair>
HilbertSpace::State::getState() const {
    return states;
}

std::string
HilbertSpace::State::getLabel() const {
    if (label_override != "") return label_override;
    std::string label = std::string();

    label += "{ ";
    for (const SPStatePair & pair : states) {
        label += pair.getLabel() + " ";
    }
    label += "}";

    return label;
}

State 
HilbertSpace::State::normalize() const {
    State newState = *this;
    double norm = (newState * newState).real();
    return newState * (1. / sqrt(norm));
}

State 
HilbertSpace::State::antisym() const {
    std::vector<SingleParticleStatePair> newStates;

    for (const SPStatePair & pair : states){
        newStates.push_back(pair);
        newStates.push_back(
            SPStatePair(
                pair.getSecondField(),
                pair.getFirstField(),
                pair.getCoef() * (-1.0)));
    }

    return State(newStates).normalize();
}

State 
HilbertSpace::State::sym() const {
    std::vector<SingleParticleStatePair> newStates;

    for (const SPStatePair & pair : states){
        newStates.push_back(pair);
        newStates.push_back(
            SPStatePair(
                pair.getSecondField(),
                pair.getFirstField(),
                pair.getCoef()));
    }

    return State(newStates).normalize();
}

State
HilbertSpace::State::operator+(const State &state) const {
    std::vector<SPStatePair> newStates;
    for (SPStatePair & pair : this->getState())
        newStates.push_back(pair);
    for (SPStatePair & pair : state.getState())
        newStates.push_back(pair);

    return State(newStates);
}

State
HilbertSpace::State::operator-(const State &state) const {
    std::vector<SPStatePair> newStates;
    for (SPStatePair & pair : this->getState())
        newStates.push_back(pair);
    for (SPStatePair & pair : state.getState())
        newStates.push_back(pair * (-1.0));

    return State(newStates);
}

State
HilbertSpace::State::operator*(Complex c) const {
    std::vector<SPStatePair> newStates;
    for (SPStatePair & pair : this->getState())
        newStates.push_back(pair * c);

    return State(newStates);
}

ComplexHighRes
HilbertSpace::State::operator*(const State &state) const {
    std::vector<SPStatePair> states1 = this->getState();
    std::vector<SPStatePair> states2 = state.getState();

    ComplexHighRes result;

    for (const SPStatePair &pair1 : states1) {
        for (const SPStatePair &pair2 : states2) {
            Complex field1Coef  = pair1.getCoef();
            SPState field1Left  = pair1.getFirstField();
            SPState field1Right = pair1.getSecondField();
            Complex field2Coef  = pair2.getCoef();
            SPState field2Left  = pair2.getFirstField();
            SPState field2Right = pair2.getSecondField();

            result += 
                (field1Left * field2Left) * 
                (field1Right * field2Right) * 
                (std::conj(field1Coef) * field2Coef);
        }
    }

    return result;
}

//////////////////////////////
//  SingleParticleOperator  //
//////////////////////////////

HilbertSpace::SingleParticleOperator::SingleParticleOperator(
    const SingleParticleFunction &left_,
    const SingleParticleFunction &right_):
    left(left_),
    right(right_) {}

State
HilbertSpace::SingleParticleOperator::operator*(
    const State &state) const {
    std::vector<SPStatePair> newStates;

    for (SPStatePair &pair : state.getState()) {
        SPState leftField  = pair.getFirstField();
        SPState rightField = pair.getSecondField();

        newStates.push_back(
            SPStatePair(
                SingleParticleState(
                    left(leftField.getField()), 
                    leftField.getSpin(),
                    leftField.getLabel()),
                SingleParticleState(
                    right(rightField.getField()), 
                    rightField.getSpin(),
                    rightField.getLabel()),
                pair.getCoef()));
    }
    return State(newStates);
}

ComplexHighRes
HilbertSpace::SingleParticleOperator::operatorValue(
    const State &left, 
    const State &right) const {
    return left * (*this * right);
}

//////////////////////////////////
// DoubleParticleScalarOperator //
//////////////////////////////////

DPSOperator::DoubleParticleScalarOperator(
    const DoubleParticleScalarFunction &func_):
    func(func_) {}

/*
State
DPSOperator::operator*(const State &) const {
    return State(std::vector());
}
*/

ComplexHighRes
DPSOperator::operatorValue(const State &left, const State &right) const {
    std::vector<SPStatePair> states1 = left.getState();
    std::vector<SPStatePair> states2 = right.getState();

    ComplexHighRes result;

    for (const SPStatePair &pair1 : states1) {
        for (const SPStatePair &pair2 : states2) {
            Complex field1Coef  = pair1.getCoef();
            SPState field1Left  = pair1.getFirstField();
            SPState field1Right = pair1.getSecondField();
            Complex field2Coef  = pair2.getCoef();
            SPState field2Left  = pair2.getFirstField();
            SPState field2Right = pair2.getSecondField();

            // TODO: It seems like Spin should be in MathUtilities, spin selection rule shouldn't be evaluated here!!!
            if  ((field1Left.getSpin() == field2Left.getSpin()) &&
                (field1Right.getSpin() == field2Right.getSpin()))
                result += twoSiteIntegral(
                    field1Left.getField(), field1Right.getField(), 
                    func, 
                    field2Left.getField(), field2Right.getField()) * 
                    (std::conj(field1Coef) * field2Coef);
        }
    }

    return result;
}

//////////////////////////////
//         Operator         //
//////////////////////////////

HilbertSpace::Operator::Operator(
    const std::vector<SingleOperator *> &operators_):
    operators(operators_) {}

HilbertSpace::Operator::~Operator() {

    //TODO: HACK!
    /*for (auto op : operators)
        delete (op);*/
}

Operator
HilbertSpace::Operator::operator+(const Operator &ops) const {
    std::vector<SingleOperator *> v1 = this->getOperator();
    std::vector<SingleOperator *> v2 = ops.getOperator();

    v1.insert(v1.end(), v2.begin(), v2.end());

    return Operator(v1);
}

std::vector<SingleOperator *>
HilbertSpace::Operator::getOperator() const {
    return operators;
}