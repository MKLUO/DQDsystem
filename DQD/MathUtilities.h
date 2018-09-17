#pragma once

#include <complex>
typedef std::complex<double> Complex;

namespace MU
{
    class ScalarField
    {
    public:
        Complex operator*(const ScalarField&) const;
    };      
}