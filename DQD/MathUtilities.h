#pragma once

#include <functional>
#include <vector>
#include <complex>

typedef std::complex<double> Complex;

typedef std::function<Complex(int, int)> SingleParticleScalarFunction;
typedef std::function<Complex(int, int, int, int)> DoubleParticleScalarFunction;

class ScalarField {
public:
    ScalarField(int, int, double, SingleParticleScalarFunction);

    ScalarField operator+(ScalarField) const;

    ScalarField operator*(Complex) const;

    Complex operator*(ScalarField) const;

    ScalarField operator*(SingleParticleScalarFunction) const;

private:
    std::vector<Complex> data;
    int width, height, gridSize;
};

typedef std::function<ScalarField(ScalarField)> SingleParticleFunction;

Complex
twoSiteIntegral(const ScalarField &, const ScalarField &, const DoubleParticleScalarFunction &, const ScalarField &,
                const ScalarField &);
