#pragma once

// TODO: Abstraction of Coordinate (x & y)

#include <functional>
#include <string>
#include <vector>
#include <complex>

constexpr double M_PI = 3.1415926535897932384626433832795028841971;

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

class ScalarField;

using SingleParticleScalarFunction = std::function<Complex(double, double)>;
using DoubleParticleScalarFunction = std::function<Complex(double, double, double, double)>;

using SingleParticleFunction = std::function<ScalarField(ScalarField)>;

using Matrix = std::vector<std::vector<ComplexHighRes>>;

struct SystemScale {
	
	SystemScale(int, int, double);

	static SystemScale defaultScale();

	const int width, height;
	const double gridSize;
};

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

	Complex value() const;

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

	void reserve(int);

private:
	std::vector<Complex> data;
};

Complex sqrt(const ComplexContainer &);
ComplexContainer operator*(const double, const ComplexContainer &);
std::ostream & operator<<(std::ostream &, const ComplexContainer &);

class Spin
{
public:
	enum class Type {
		Up,
		Down,
		None
	};

	explicit
	Spin(Type);

	double
	operator*(const Spin &)
	const;

	bool 
	operator==(const Spin &)
	const;

	Type
	getType()
	const;

	std::string
	getLabel()
	const;

private:
	const Type type;
};

class ScalarField {
public:

	explicit
	ScalarField(
		const SystemScale);

	explicit
	ScalarField(
		const SystemScale, 
		const std::vector<Complex> &);

	explicit
	ScalarField(
		const SystemScale, 
		const SingleParticleScalarFunction &);

	// Operators

	ScalarField operator+(const ScalarField &) const;

	ScalarField operator-(const ScalarField &) const;

	ScalarField operator*(const SingleParticleScalarFunction &) const;

	ScalarField operator*(const SingleParticleFunction &) const;

	ScalarField operator*(Complex) const;

	ScalarField operator/(Complex) const;

	ScalarField operator*(double) const;

	ComplexHighRes operator*(const ScalarField &) const;

	ScalarField operator^(const ScalarField &) const;

	ScalarField conj() const;

	// Access

	std::vector<Complex> getDatas() const;

	Complex getData(int, int) const;

	double getX(int) const;

	double getY(int) const;

	int getIndex(int, int) const;

	SystemScale getScale() const;

	std::vector<double> norm() const;

	// Modify
	
	Complex &Data(int, int);

	// Utils

	ScalarField normalize() const;

	void plotTemp() const; // For DEBUG

private:
	const SystemScale scale;
	std::vector<Complex> data;
};

ScalarField operator*(Complex, const ScalarField &);

ScalarField operator*(double, const ScalarField &);

ScalarField operator*(const SingleParticleScalarFunction &, const ScalarField &);

ScalarField operator*(const SingleParticleFunction &, const ScalarField &);

// Math Utilities

std::vector<Complex>
functionToField(
	const SingleParticleScalarFunction &, 
	const SystemScale);

int
fieldIndex(
	const int i,
	const int j,
	const int width,
	const int height
);

double
fieldCoord(
	const int i,
	const int width,
	const double gridSize
);

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




