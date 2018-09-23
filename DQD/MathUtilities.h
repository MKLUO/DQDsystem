#pragma once

#include <functional>
#include <vector>
#include <complex>
typedef std::complex<double> Complex;

class ScalarField;

typedef std::function<Complex(double, double)> SingleParticleScalarFunction;
typedef std::function<Complex(double, double, double, double)> DoubleParticleScalarFunction;

typedef std::function<ScalarField(ScalarField)> SingleParticleFunction;

class ScalarField {
public:
    ScalarField(int, int, double, const std::vector<Complex> &);

    ScalarField(int, int, double, const SingleParticleScalarFunction &);

    // Operators

    ScalarField operator+(const ScalarField &) const;

    ScalarField operator*(Complex) const;

    Complex operator*(const ScalarField &) const;

    ScalarField operator*(const SingleParticleScalarFunction &) const;

    ScalarField operator*(const SingleParticleFunction &) const;

    // Access

    Complex at(int, int) const;

    double getX(int) const;

    double getY(int) const;

    int getIndex(int, int) const;

    int getWidth() const;

    int getHeight() const;

    double getGridSize() const;

    // Math Utilities

    ScalarField laplacian() const;

    ScalarField angularMomentum() const;

private:
    std::vector<Complex> data;
    int width, height;
    double gridSize;
};

Complex
twoSiteIntegral(const ScalarField &, const ScalarField &, const DoubleParticleScalarFunction &, const ScalarField &,
                const ScalarField &);
