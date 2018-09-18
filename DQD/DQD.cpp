#include "DQD.h"

// Abbreviations
typedef DQD::HilbertSpace::Hamiltonian                      Hamiltonian;
typedef DQD::HilbertSpace::State                            State;
typedef DQD::HilbertSpace::State::TwoParticleWaveFunction   Field;
typedef std::vector<Field>                                  Fields;

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
    MU::ScalarField field11 = field1.getField1();
    MU::ScalarField field12 = field1.getField2();
    MU::ScalarField field21 = field2.getField1();
    MU::ScalarField field22 = field2.getField2();

    Complex result;

    // Kinetic energy
    result += ham.kineticConstant() * field11 * MU::laplacian(field21);
    result += ham.kineticConstant() * field12 * MU::laplacian(field22);

    // External Potential
    result += ham.coulombConstant() * field11 * MU::inverseR(field21);
    result += ham.coulombConstant() * field12 * MU::inverseR(field22);

    // Internal Coulomb energy
    result += ham.coulombConstant() * MU::operate(field11, field12, MU::inverseR(), field21, field22);

    return result;
}

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





