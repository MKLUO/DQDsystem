#pragma once

#include <functional>

typedef std::complex<double> Complex;

typedef std::function<Complex(int, int)>           SingleParticleScalarFunction;
typedef std::function<Complex(int, int, int, int)> DoubleParticleScalarFunction;

class ScalarField {
public:
    ScalarField(int, int);

    void setDataWithFunction(Complex (*)(int, int));

    // Operators
    // "*": Inner product of two ScalarField
    // "^": Tensor product of two ScalarField
    Complex operator*(const ScalarField&) const;

private:
    std::vector<Complex> data;
    int width;
};

typedef std::function<ScalarField(ScalarField)>    SingleParticleFunction;

ScalarField scalarField();

ScalarField laplacian(const ScalarField&);
ScalarField multiply(const ScalarField&, const ScalarField&);

// Note: This should be generalized
Complex twoSiteInverseRIntegral(const ScalarField&, const ScalarField&, const ScalarField&, const ScalarField&);
