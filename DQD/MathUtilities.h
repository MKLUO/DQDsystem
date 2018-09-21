#pragma once

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

class ScalarFunction {
public:
    ScalarFunction(Complex (*)(int, int));
};

class SingleParticleFunction {

};

class DoubleParticleFunction {

};

ScalarField scalarField();

ScalarField laplacian(const ScalarField&);
ScalarField multiply(const ScalarField&, const ScalarField&);

// Note: This should be generalized
Complex twoSiteInverseRIntegral(const ScalarField&, const ScalarField&, const ScalarField&, const ScalarField&);

Complex fockDarwin(int, int);