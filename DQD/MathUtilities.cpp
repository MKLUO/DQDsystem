#include "DQD.h"

// Abbreviations
typedef DQD::HilbertSpace::MathUtilities					MU;
typedef DQD::HilbertSpace::MathUtilities::ScalarField		ScalarField;

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
	return ScalarField(width * height);
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

ScalarField::ScalarField(int size)
{
	data.resize(size);
}

Complex ScalarField::operator*(const ScalarField&) const
{
    Complex result;

    // TODO: integration (Single particle wavefunction inner product);

    return result;
}