#pragma once

#include <vector>

#include "MathUtilities.h"

class DQD
{
public:
    class HilbertSpace
    {
    public:
        class State;
        class Hamiltonian;

        class Hamiltonian
        {
        public:
            // Operators
            // "*": Hamiltonian operate on State
            State operator*(const State&) const;
            //State::TwoParticleWaveFunction operator*(const State::TwoParticleWaveFunction&) const;

        private:
            // Parameters of terms in Hamiltonian (Kinetic, potential, repulsion)
        };

        class State
        {
        public:
            class TwoParticleWaveFunction
            {
            public:
                MU::ScalarField getField1() const;
                MU::ScalarField getField2() const;

                // Operators
                // "*": Inner product of two TwoParticleWaveFunctions
                Complex operator*(const TwoParticleWaveFunction&) const;

            private:
                MU::ScalarField scalarField1;
                MU::ScalarField scalarField2;
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

    private:
        int width, height;
    };

    DQD();
    ~DQD();

private:
};