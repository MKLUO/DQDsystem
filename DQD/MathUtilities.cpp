#include "MathUtilities.h"

using MU::ScalarField;

//////////////////////////////
//      	MU        		//
//////////////////////////////

ScalarField laplacian(const ScalarField& field)
{
	ScalarField newField;
	
	return newField;
}

ScalarField inverseR(const ScalarField&)
{
	ScalarField newField;
	
	return newField;
}

Complex twoSiteInverseRIntegral(const ScalarField&, const ScalarField&, const ScalarField&, const ScalarField&)
{
	Complex result;
	
	return result;
}


//////////////////////////////
//      ScalarField        	//
//////////////////////////////

Complex ScalarField::operator*(const ScalarField&) const
{
    Complex result;

    // TODO: integration (Single particle wavefunction inner product);

    return result;
}