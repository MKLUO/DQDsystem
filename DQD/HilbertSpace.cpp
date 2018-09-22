#include <utility>

#include "HilbertSpace.h"

// Abbreviations
typedef SingleParticleScalarFunction                SPSFunction;
typedef DoubleParticleScalarFunction                DPSFunction;

typedef SingleParticleFunction                      SPFunction;

typedef HilbertSpace::SingleParticleState           SPState;
typedef HilbertSpace::SingleParticleStatePair       SPStatePair;
typedef HilbertSpace::State                         State;

typedef HilbertSpace::SingleOperator                SingleOperator;
typedef HilbertSpace::SingleParticleOperator        SPOperator;
typedef HilbertSpace::DoubleParticleScalarOperator  DPSOperator;
typedef HilbertSpace::Operator                      Operator;


//////////////////////////////
//       HilbertSpace       //
//////////////////////////////

HilbertSpace::HilbertSpace(int width_, int height_) {
    width   = width_;
    height  = height_;
}

SPState HilbertSpace::createSingleParticleState(const SPSFunction& function) const {

}
Operator HilbertSpace::createOperator(const SingleParticleFunction& function1, const SingleParticleFunction& function2) const {

}
Operator HilbertSpace::createOperator(const DoubleParticleScalarFunction& function) const {

}

Complex HilbertSpace::operatorValue(const State& stateLeft, const Operator& ops, const State& stateRight) const {

    Complex result;

    for (SingleOperator* op : ops.getOperator()) {
        result += op->operatorValue(stateLeft, stateRight);
    }

    return result;
}

//////////////////////////////
//          State           //
//////////////////////////////

State::State()
{

}

State::State(Fields fields_)
{
    fields = std::move(fields_);
}

Fields State::getFields() const
{
    return fields;
}

State State::operator+(const State& state) const
{
    Fields fields1 = this->getFields();
    Fields fields2 = state.getFields();

    fields1.insert(fields1.end(), fields2.begin(), fields2.end());

    return State(fields1);
}

Complex State::operator*(const State& state) const
{
    Fields fields1 = this->getFields();
    Fields fields2 = state.getFields();

    Complex result;

    for (Field field1 : fields1)
        for (Field field2 : fields2)
            result += field1 * field2;

    return result;
}

//////////////////////////////
// TwoParticleWaveFunction  // * Abbreviated as "Field"
//////////////////////////////

ScalarField Field::field1() const
{
    return scalarField1;
}

ScalarField Field::field2() const
{
    return scalarField2;
}

Complex Field::operator*(const Field& func) const
{
    return this->field1() * func.field1() +
            this->field2() * func.field2();
}





