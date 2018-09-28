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

    ScalarField operator-(const ScalarField &) const;

    ScalarField operator*(Complex) const;

    Complex operator*(const ScalarField &) const;

    ScalarField operator*(const SingleParticleScalarFunction &) const;

    ScalarField operator*(const SingleParticleFunction &) const;

    // Access

    std::vector<Complex> getDatas() const;

    Complex getData(int, int) const;

    double getX(int) const;

    double getY(int) const;

    int getIndex(int, int) const;

    int getWidth() const;

    int getHeight() const;

    double getGridSize() const;


private:
    Complex& at(int, int);

    std::vector<Complex> data;
    int width, height;
    double gridSize;
};

// Math Utilities

Complex
twoSiteIntegral(const ScalarField &, const ScalarField &, const DoubleParticleScalarFunction &, const ScalarField &,
                const ScalarField &);

ScalarField
laplacian(const ScalarField &);

ScalarField
angularMomentum(const ScalarField &);

extern SingleParticleScalarFunction
x_func;

extern SingleParticleScalarFunction
y_func;

Complex
gaussian(double, double, double);


