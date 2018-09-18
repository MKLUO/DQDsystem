#ifndef MATHUTILITIES_H
#define MATHUTILITIES_H

#include <complex>
typedef std::complex<double> Complex;

namespace MU
{
    class ScalarField
    {
    public:
		Complex operator*(const ScalarField&) const;
    };  

	ScalarField laplacian(const ScalarField&);
	ScalarField inverseR(const ScalarField&);
	
	// Note: This should be generalized
	Complex twoSiteInverseRIntegral(const ScalarField&, const ScalarField&, const ScalarField&, const ScalarField&);
}

#endif //MATHUTILITIES_H