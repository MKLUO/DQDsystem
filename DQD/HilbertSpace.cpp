#include "HilbertSpace.h"
#include <utility>

// Abbreviations
typedef HilbertSpace::SingleParticleState       SPState;
typedef HilbertSpace::SingleParticleStatePair   SPStatePair;
typedef HilbertSpace::State                     State;

typedef HilbertSpace::SingleOperator            SingleOperator;
typedef HilbertSpace::SingleParticleOperator    SPOperator;
typedef HilbertSpace::DoubleParticleOperator    DPOperator;
typedef HilbertSpace::Operator                  Operator;


//////////////////////////////
//       HilbertSpace       //
//////////////////////////////

HilbertSpace::HilbertSpace(int width_, int height_)
{
    width = width_;
    height = height_;
}

Complex HilbertSpace::product(const State& state1, const Hamiltonian& ham, const State& state2) const
{
    Fields fields1 = state1.getFields();
    Fields fields2 = state2.getFields();

    Complex result;

    for (Field field1 : fields1)
        for (Field field2 : fields2)
            result += product(field1, ham, field2);

    return result;
}

Complex HilbertSpace::product(const Field& field1, const Hamiltonian& ham, const Field& field2) const
{
    ScalarField field11 = field1.field1();
    ScalarField field12 = field1.field2();
    ScalarField field21 = field2.field1();
    ScalarField field22 = field2.field2();

    Complex result;

    // Kinetic energy
    result += ham.kineticConstant() * (field11 * mu->laplacian(field21));
    result += ham.kineticConstant() * (field12 * mu->laplacian(field22));

    // External Potential
    result += field11 * mu->multiply(ham.getPotential(), field21);
    result += field12 * mu->multiply(ham.getPotential(), field22);

    // Internal Coulomb energy
    result += ham.coulombConstant() * mu->twoSiteInverseRIntegral(field11, field12, field21, field22);

    return result;
}

//////////////////////////////
//      Hamiltonian         //
//////////////////////////////

Complex Hamiltonian::kineticConstant() const
{
    return Complex();
}

Complex Hamiltonian::coulombConstant() const
{
    return Complex();
}

ScalarField Hamiltonian::getPotential() const
{
    return potential;
}

State Hamiltonian::operator*(const State& state) const
{
    Fields fields = state.getFields();

    for (Field& field : fields)
        //field = (*this) * field;

    return State(fields);
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





