#pragma once

#include <vector>

#include <complex>
typedef std::complex<double> Complex;

class DQD
{
public:
    struct Setting
    {
        static Setting defaultSetting();

        int width, height;
        double omega0;
    };

    class HilbertSpace
    {
    public:
        class MathUtilities;
        class Hamiltonian;
        class State;

        class MathUtilities
        {
        public:
            class ScalarField
            {
            public:
                ScalarField(int);
                Complex operator*(const ScalarField&) const;
            private:
                std::vector<Complex> data;
            };

            MathUtilities(int, int);

            ScalarField scalarField() const;

            ScalarField laplacian(const ScalarField&) const;
            ScalarField multiply(const ScalarField&, const ScalarField&) const;

            // Note: This should be generalized
            Complex twoSiteInverseRIntegral(const ScalarField&, const ScalarField&, const ScalarField&, const ScalarField&) const;

        private:
            int width, height;
        };

        class Hamiltonian
        {
        public:
            Complex kineticConstant() const;
			Complex coulombConstant() const;

            MathUtilities::ScalarField getPotential() const;
			
			// Operators
			// "*": Operation of Hamiltonian on State
			State operator*(const State&) const;
				
        private:
            MathUtilities::ScalarField potential;
            // Parameters of terms in Hamiltonian (Kinetic, potential, repulsion)
        };

        class State
        {
        public:
            class TwoParticleWaveFunction
            {
            public:
                MathUtilities::ScalarField getField1() const;
                MathUtilities::ScalarField getField2() const;

                // Operators
                // "*": Inner product of two TwoParticleWaveFunctions
                Complex operator*(const TwoParticleWaveFunction&) const;

            private:
                MathUtilities::ScalarField scalarField1;
                MathUtilities::ScalarField scalarField2;
            };

            State();
            State(std::vector<TwoParticleWaveFunction>);

            std::vector<TwoParticleWaveFunction> getFields() const;

            // Operators
            // "+": Addition of States
            // "*": Inner product of two States
            State operator+(const State&) const;
            Complex operator*(const State&) const;

        private:
            std::vector<TwoParticleWaveFunction> fields;
        };

        HilbertSpace(Setting);

        Complex product(const State&, const Hamiltonian&, const State&) const;
        Complex product(const State::TwoParticleWaveFunction&, const Hamiltonian&, const State::TwoParticleWaveFunction&) const;

    private:
        MathUtilities* mu;
        Setting setting;
    };

    DQD();
    ~DQD();

private:
    HilbertSpace* hilbertSpace;
};