#pragma once

// TODO: Abstraction of Coordinate (x & y)

#include <functional>
#include <string>
#include <vector>
#include <complex>

// For i
using namespace std::literals::complex_literals; 

using Complex = std::complex<double>;
class ComplexContainer;
const int COMPLEX_MAX_SIZE = 50000; // TODO: What is a reasonable max size?
const double COMPLEX_SHRINK_RATIO = 0.8;

// ========================================================
// ComplexHighRes can chosen as ComplexContainer or Complex
// ========================================================
#define HIGH_RES_COMPLEX

#ifdef HIGH_RES_COMPLEX
	using ComplexHighRes = ComplexContainer;
#else
	using ComplexHighRes = Complex;
#endif

enum class Spin {
	Up,
	Down,
	None
};

std::string spinSign(const Spin &);

class ScalarField;

using SingleParticleScalarFunction = std::function<Complex(double, double)>;
using DoubleParticleScalarFunction = std::function<Complex(double, double, double, double)>;

using SingleParticleFunction = std::function<ScalarField(ScalarField)>;

using Matrix = std::vector<std::vector<ComplexHighRes>>;

namespace PhysicsContant {

	// TODO: Physics related entities shouldn't be here

	const double e = 1.60217662e-19;
	const double me = 9.10938356e-31;
	const double hBar = 1.0545718e-34;
	const double epsilon = 8.854187817e-12;

	// For Silicon
	const double m = me * 0.191;
	const double kappa = 7.64;
}

// This ComplexContainer class represents a complex number. While it performs addition or substraction, instead of directly evaluating the result, it stores a list of std::complex<double> values according to the input values. "data" may grow big so a check on every operation is a must.
class ComplexContainer {
public:

	ComplexContainer();

	ComplexContainer(Complex);

	ComplexContainer(double, double);

	ComplexContainer(double);

	ComplexContainer(std::vector<Complex>);

	// operator Complex() const;

	ComplexContainer operator+(const ComplexContainer &) const;
	ComplexContainer operator-(const ComplexContainer &) const;
	ComplexContainer operator*(const ComplexContainer &) const;
	ComplexContainer operator/(const ComplexContainer &) const;

	ComplexContainer operator-() const;

	double real() const;

	double imag() const;

	double norm() const;

	ComplexContainer conj() const;

	int size() const;

	bool isZero() const;
	
	// Non-const methods.

	void operator+=(const ComplexContainer &);
	void operator+=(const Complex &);

	void shrink(double);

	void shrink(int);

	// Utils

	void reserve(int);

private:
	// Debug only
	Complex value() const;

	std::vector<Complex> data;
};

ComplexContainer operator*(double, const ComplexContainer &);
std::ostream & operator<<(std::ostream &, const ComplexContainer &);


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

	ComplexHighRes operator*(const ScalarField &) const;

	ScalarField operator^(const ScalarField &) const;

	ScalarField conj() const;

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

ComplexHighRes
twoSiteIntegral(const ScalarField &, const ScalarField &, const DoubleParticleScalarFunction &, const ScalarField &,
				const ScalarField &);

ScalarField
reverse(const ScalarField&);

double 
oneMinus_sqrtOneMinusXX_divideX(const double&);

// Complex
// value(const ComplexHighRes&);

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
planeWave(double, double);

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




