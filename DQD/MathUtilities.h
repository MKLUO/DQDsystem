#pragma once

// TODO: Abstraction of Coordinate (x & y)

// TODO: design namespace

#include <functional>
#include <vector>
#include <complex>

using Complex = std::complex<double>;
using namespace std::complex_literals;

class ScalarField;

using SingleParticleScalarFunction = std::function<Complex(double, double)>;
using DoubleParticleScalarFunction = std::function<Complex(double, double, double, double)>;

using SingleParticleFunction = std::function<ScalarField(ScalarField)>;

namespace Physics {
    // TODO
    const double e = 1.60217662e-19;
    const double m = 9.10938356e-31 * 0.191;
    const double hBar = 1.0545718e-34;
    const double epsilon = 8.854187817e-12;
}

class ScalarField {
public:
    ScalarField(int, int, double);

    ScalarField(int, int, double, const std::vector<Complex> &);

    ScalarField(int, int, double, const SingleParticleScalarFunction &);

    // Operators

    ScalarField operator+(const ScalarField &) const;

    ScalarField operator-(const ScalarField &) const;

    ScalarField operator*(const SingleParticleScalarFunction &) const;

    ScalarField operator*(const SingleParticleFunction &) const;

    ScalarField operator*(Complex) const;

    ScalarField operator*(double) const;

    Complex operator*(const ScalarField &) const;

    ScalarField operator^(const ScalarField &) const;

    // Access

    std::vector<Complex> getDatas() const;

    Complex getData(int, int) const;

    Complex &setData(int, int);

    double getX(int) const;

    double getY(int) const;

    int getIndex(int, int) const;

    int getWidth() const;

    int getHeight() const;

    double getGridSize() const;

    std::vector<double> norm() const;

private:


    std::vector<Complex> data;
    int width, height;
    double gridSize;
};

ScalarField operator*(Complex, const ScalarField &);

ScalarField operator*(double, const ScalarField &);

ScalarField operator*(const SingleParticleScalarFunction &, const ScalarField &);

ScalarField operator*(const SingleParticleFunction &, const ScalarField &);

// Math Utilities

Complex
twoSiteIntegral(const ScalarField &, const ScalarField &, const DoubleParticleScalarFunction &, const ScalarField &,
                const ScalarField &);

ScalarField
reverse(const ScalarField&);

// ScalarFields

extern SingleParticleScalarFunction
        x_field;

extern SingleParticleScalarFunction
        y_field;

extern SingleParticleScalarFunction
        xx_field;

extern SingleParticleScalarFunction
        yy_field;

extern SingleParticleScalarFunction
        sho_field;

// ScalarFields with settings needed

SingleParticleScalarFunction
scalar(Complex);

SingleParticleScalarFunction
gaussian(double);

SingleParticleScalarFunction
quartic(double);

DoubleParticleScalarFunction
rInv_field(double);

// ScalarFunctions

extern SingleParticleFunction
        identity;

extern SingleParticleFunction
        laplacian;

extern SingleParticleFunction
        angularMomentum;




