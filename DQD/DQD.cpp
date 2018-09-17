#include "DQD.h"

// Abbreviations
typedef DQD::HilbertSpace::Hamiltonian Hamiltonian;
typedef DQD::HilbertSpace::State State;
typedef DQD::HilbertSpace::State::TwoParticleWaveFunction Field;
typedef std::vector<Field> Fields;

//////////////////////////////
//          DQD             //
//////////////////////////////

DQD::DQD()
{
}


DQD::~DQD()
{
}

//////////////////////////////
//      HilbertSpace        //
//////////////////////////////

//////////////////////////////
//      Hamiltonian         //
//////////////////////////////

State Hamiltonian::operator*(const State& state) const
{
    Fields fields = state.getFields();

    for (Field& field : fields)
        //field = (*this) * field;

    return State(fields);
}
/*
Field Hamiltonian::operator*(const Field& field) const
{
    
}
*/
//////////////////////////////
//          State           //
//////////////////////////////

State::State()
{
    
}

State::State(Fields fields_)
{
    fields = fields_;
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

MU::ScalarField Field::getField1() const
{
    return scalarField1;
}

MU::ScalarField Field::getField2() const
{
    return scalarField2;
}

Complex Field::operator*(const Field& func) const
{
    return  this->getField1() * func.getField1() +
            this->getField2() * func.getField2();
}





