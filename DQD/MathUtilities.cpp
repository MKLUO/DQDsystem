#include "HilbertSpace.h"

// Abbreviations
typedef HilbertSpace::MathUtilities					MU;
typedef HilbertSpace::MathUtilities::ScalarField		ScalarField;

typedef HilbertSpace::State                            State;
typedef HilbertSpace::State::TwoParticleWaveFunction   Field;
typedef std::vector<Field>                                Fields;

//////////////////////////////
//      	MU        		//
//////////////////////////////

MU::MathUtilities(int width_, int height_)
{
	width = width_;
	height = height_;
}

ScalarField MU::scalarField() const
{
	return ScalarField(width, height);
}

ScalarField MU::laplacian(const ScalarField& field) const
{
	ScalarField newField = scalarField();
	
	return newField;
}

ScalarField MU::multiply(const ScalarField& func, const ScalarField& field) const
{
	ScalarField newField = scalarField();
	
	return newField;
}

Complex MU::twoSiteInverseRIntegral(const ScalarField&, const ScalarField&, const ScalarField&, const ScalarField&) const
{
	Complex result;
	
	return result;
}


//////////////////////////////
//      ScalarField        	//
//////////////////////////////

ScalarField::ScalarField(int width_, int height_)
{
	width = width_;
	data.resize(width_ * height_);
}

void ScalarField::setDataWithFunction(Complex (*func)(int, int))
{
	for (int i = 0; i < data.size(); i++) {
		int x = i % width;
		int y = (i - x) / width;
		data[i] = func(x, y);
	}
}

Complex ScalarField::operator*(const ScalarField&) const
{
    Complex result;

    // TODO: integration (Single particle wavefunction inner product);

    return result;
}

State ScalarField::operator^(const ScalarField& newField) const
{
	Field field = Field(*this, newField);
	return State(Fields{field});
}

//////////////////////////////
//        Utilities        	//
//////////////////////////////

Complex MU::fockDarwin(int x, int y)
{
	return Complex();
}
