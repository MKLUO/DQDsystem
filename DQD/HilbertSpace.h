#pragma once

#include <vector>
#include <complex>

#include "MathUtilities.h"


// HilbertSpace: Implements 2D two particle Hilbert space with finite size.
class HilbertSpace {
public:
    // Pre-declaration of classes
    // Structure:
    //      State
    //          `- SingleParticleStatePair
    //                  `- SingleParticleState
    //
    //      Operator
    //          `- (SingleOperator*) SingleParticleOperator
    //              `- SingleParticleFunction
    //          `- (SingleOperator*) DoubleParticleScalarOperator
    //              `- DoubleParticleScalarFunction
    //
    class   SingleParticleState;
    class   SingleParticleStatePair;
    class   State;

    class   SingleOperator;
    class   SingleParticleOperator;
    class   DoubleParticleOperator;
    class   Operator;

    // SingleParticleState: Represents a one-particle wavefunction.
    class SingleParticleState {
    public:
        SingleParticleState(const ScalarField&);

        ScalarField getField() const;

        // Operators
        // "+": Addition of two SingleParticleState
        // "*": Inner product of two SingleParticleState
        // "^": Tensor product of two SingleParticleState
        SingleParticleState operator+(const SingleParticleState&) const;
        Complex operator*(const SingleParticleState&) const;
        SingleParticleStatePair operator^(const SingleParticleState&) const;

    private:
        ScalarField field;
    };

    // SingleParticleStatePair: Represents a separable two-particle wavefunction.
    class SingleParticleStatePair {
    public:
        SingleParticleStatePair(const SingleParticleState&, const SingleParticleState&);

        ScalarField getFirstField() const;
        ScalarField getSecondField() const;

    private:
        SingleParticleState first;
        SingleParticleState second;
    };

    // State: Represents a general two-particle wavefunction.
    class State {
    public:
        State(const std::vector<SingleParticleStatePair>&);
        ~State();

        std::vector<SingleParticleStatePair> getState() const;

        // Operators
        // "+": Addition of two State
        // "*": Multiply with a scalar
        // "*": Inner product of two State
        State operator+(const State&) const;
        State operator*(Complex) const;
        Complex operator*(const State&) const;

        // TODO: print()

    private:
        std::vector<SingleParticleStatePair> states;
    };

    // Operator: A interface for both Single/DoubleParticlOoperator.
    class SingleOperator {
    public:
        virtual State operator*(const State&) const = 0;

        virtual Complex operatorValue(const State&, const State&) const = 0;
    };

    // SingleParticleOperator: Represents a one-particle operator.
    class SingleParticleOperator: public SingleOperator {
    public:
        SingleParticleOperator(const SingleParticleFunction&, const SingleParticleFunction&);

        // Operators
        // "*": Operation of SingleParticleOperator on State
        State operator*(const State&) const override;

        Complex operatorValue(const State&, const State&) const override;

    private:
        SingleParticleFunction left;
        SingleParticleFunction right;
    };

    // DoubleParticleScalarOperator: Represents a two-particle scalar operator.
    class DoubleParticleScalarOperator: public SingleOperator {
    public:
        DoubleParticleScalarOperator(const DoubleParticleScalarFunction&);

        // Operators
        // "*": Operation of DoubleParticleOperator on State
        State operator*(const State&) const override;

        Complex operatorValue(const State&, const State&) const override;

    private:
        DoubleParticleScalarFunction func;
    };

    class Operator {
    public:
        Operator(const std::vector<SingleOperator*>&);
        ~Operator();

        // Operators
        // "+": Composition of two Operators
        // "*": Operation of Operator on State
        Operator operator+(const Operator &) const;
        //State operator*(const State &) const;

        std::vector<SingleOperator*> getOperator() const;

        // TODO: How to handle generic operations?

    private:
        std::vector<SingleOperator*> operators;
    };

    // Constructor input: width and height of the Hilbert space.
    HilbertSpace(int, int);

    SingleParticleState createSingleParticleState(const SingleParticleScalarFunction&) const;
    Operator createOperator(const SingleParticleFunction&, const SingleParticleFunction&) const;
    Operator createOperator(const DoubleParticleScalarFunction&) const;

    // TODO: proper naming of this function.
    Complex operatorValue(const State&, const Operator&, const State&) const;
    double expectationValue(const State&, const Operator&) const;

private:
    int width, height;
};