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
    // NOTE:

    class SingleParticleState;

    class State;

    class SingleOperator;

    class SingleParticleOperator;

    class DoubleParticleScalarOperator;

    class Operator;

    // SingleParticleState: Represents a one-particle wavefunction.
    class SingleParticleState {
    public:
        explicit SingleParticleState(const ScalarField &);

        ScalarField getField() const;

        // Operators
        // "+": Addition of two SingleParticleState
        // "*": Multiply with a scalar
        // "*": Inner product of two SingleParticleState
        // "^": Tensor product of two SingleParticleState
        SingleParticleState operator+(const SingleParticleState &) const;

        SingleParticleState operator*(Complex) const;

        Complex operator*(const SingleParticleState &) const;

        State operator^(const SingleParticleState &) const;

    private:
        ScalarField field;
    };

    // SingleParticleStatePair: Represents a separable two-particle wavefunction.
    class SingleParticleStatePair {
    public:
        SingleParticleStatePair(const SingleParticleState &, const SingleParticleState &);

        SingleParticleState getFirstField() const;

        SingleParticleState getSecondField() const;

        // Operators
        // "*": Multiply with a scalar
        SingleParticleStatePair operator*(Complex) const;

    private:
        SingleParticleState first;
        SingleParticleState second;
    };

    // State: Represents a general two-particle wavefunction.
    class State {
    public:

        explicit State(const SingleParticleState &, const SingleParticleState &);

        explicit State(const std::vector<SingleParticleStatePair> &);

        std::vector<SingleParticleStatePair> getState() const;

        State normalize() const;

        // Operators
        // "+": Addition of two State
        // "*": Multiply with a scalar
        // "*": Inner product of two State
        State operator+(const State &) const;

        State operator-(const State &) const;

        State operator*(Complex) const;

        Complex operator*(const State &) const;

        // TODO: print(), normalization

    private:
        std::vector<SingleParticleStatePair> states;
    };

    // Operator: A interface for both Single/DoubleParticlOoperator.
    class SingleOperator {
    public:
        //virtual State operator*(const State &) const = 0;

        virtual Complex operatorValue(const State &, const State &) const = 0;
    };

    // SingleParticleOperator: Represents a one-particle operator.
    class SingleParticleOperator : public SingleOperator {
    public:
        SingleParticleOperator(const SingleParticleFunction &, const SingleParticleFunction &);

        // Operators
        // "*": Operation of SingleParticleOperator on State
        State operator*(const State &) const;

        SingleParticleStatePair operator*(const SingleParticleStatePair &) const;

        Complex operatorValue(const State &, const State &) const override;

    private:
        SingleParticleFunction left;
        SingleParticleFunction right;
    };

    // DoubleParticleScalarOperator: Represents a two-particle scalar operator.
    class DoubleParticleScalarOperator : public SingleOperator {
    public:
        explicit DoubleParticleScalarOperator(const DoubleParticleScalarFunction &);

        // Operators
        // "*": Operation of DoubleParticleScalarOperator on State
        //State operator*(const State &) const override;

        Complex operatorValue(const State &, const State &) const override;

    private:
        DoubleParticleScalarFunction func;
    };

    class Operator {
    public:
        explicit Operator(const std::vector<SingleOperator *> &);

        ~Operator();

        // Operators
        // "+": Composition of two Operators
        // "*": Operation of Operator on State
        Operator operator+(const Operator &) const;
        //State operator*(const State &) const;

        std::vector<SingleOperator *> getOperator() const;

        // TODO: How to handle generic operations?

    private:
        std::vector<SingleOperator *> operators;
    };

    // Constructor input: width and height of the Hilbert space.
    HilbertSpace(int, int, double);

    SingleParticleState createSingleParticleState(const SingleParticleScalarFunction &) const;

    ScalarField createScalarField() const;

    ScalarField createScalarField(const std::vector<Complex> &) const;

    ScalarField createScalarField(const SingleParticleScalarFunction &) const;

    Operator createOperator(const SingleParticleFunction &, const SingleParticleFunction &) const;

    Operator createOperator(const DoubleParticleScalarFunction &) const;

    // TODO: proper naming of this function.
    Complex operatorValue(const State &, const Operator &, const State &) const;

    double expectationValue(const State &, const Operator &) const;

private:
    int width, height;
    double gridSize;
};

// Abbreviations
using SPSFunction   = SingleParticleScalarFunction;
using DPSFunction   = DoubleParticleScalarFunction;

using SPFunction    = SingleParticleFunction;

using SPState       = HilbertSpace::SingleParticleState;
using SPStatePair   = HilbertSpace::SingleParticleStatePair;
using State         = HilbertSpace::State;

using SingleOperator    = HilbertSpace::SingleOperator;
using SPOperator        = HilbertSpace::SingleParticleOperator;
using DPSOperator       = HilbertSpace::DoubleParticleScalarOperator;
using Operator          = HilbertSpace::Operator;