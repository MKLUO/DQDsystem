#pragma once

#include <string> 
#include <vector>
#include <complex>

#include "math_utilities.hpp"

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

    struct SystemScale {
        static SystemScale defaultScale();

        int width, height;
        double gridSize;
    };

    class SingleParticleState;

    class State;

    class SingleOperator;

    class SingleParticleOperator;

    class DoubleParticleScalarOperator;

    class Operator;

    // SingleParticleState: Represents a one-particle wavefunction.
    class SingleParticleState {
    public:
        explicit SingleParticleState(const ScalarField &, const Spin &, const std::string &);
        explicit SingleParticleState(const ScalarField &, const Spin &);

        ScalarField getField() const;

        Spin getSpin() const;

        std::string getLabel() const;

        // Operators
        // "+": Addition of two SingleParticleState
        // "*": Multiply with a scalar
        // "*": Inner product of two SingleParticleState
        // "^": Tensor product of two SingleParticleState
        SingleParticleState operator+(const SingleParticleState &) const;

        SingleParticleState operator-(const SingleParticleState &) const;

        SingleParticleState operator*(Complex) const;

        SingleParticleState operator/(Complex) const;

        ComplexHighRes operator*(const SingleParticleState &) const;

        SingleParticleState operator*(const SingleParticleFunction &) const;

        State operator^(const SingleParticleState &) const;

    private:
        ScalarField field;
        Spin spin;
        std::string label;

        static int auto_index;
    };

    // SingleParticleStatePair: Represents a separable two-particle wavefunction.
    class SingleParticleStatePair {
    public:
        explicit SingleParticleStatePair(const SingleParticleState &, const SingleParticleState &, const std::string & label = "");

        SingleParticleState getFirstField() const;

        SingleParticleState getSecondField() const;

        std::string getLabel() const;

        // Operators
        // "*": Multiply with a scalar
        SingleParticleStatePair operator*(Complex) const;

    private:
        SingleParticleState first;
        SingleParticleState second;
        std::string label_override;
    };

    // State: Represents a general two-particle wavefunction.
    class State {
    public:

        explicit State(const SingleParticleState &, const SingleParticleState &, const std::string & label = "");

        explicit State(const std::vector<SingleParticleStatePair> &, const std::string & label = "");

        std::vector<SingleParticleStatePair> getState() const;

        std::string getLabel() const;

        State normalize() const;

        State antisym() const;

        State sym() const;

        // Operators
        // "+": Addition of two State
        // "*": Multiply with a scalar
        // "*": Inner product of two State
        State operator+(const State &) const;

        State operator-(const State &) const;

        State operator*(Complex) const;

        ComplexHighRes operator*(const State &) const;

    private:
        std::vector<SingleParticleStatePair> states;
        std::string label_override;
    };

    // Operator: A interface for both Single/DoubleParticlOoperator.
    class SingleOperator {
    public:
        //virtual State operator*(const State &) const = 0;

        virtual ComplexHighRes operatorValue(const State &, const State &) const = 0;
    };

    // SingleParticleOperator: Represents a one-particle (separable) operator.
    class SingleParticleOperator : public SingleOperator {
    public:
        explicit SingleParticleOperator(const SingleParticleFunction &, const SingleParticleFunction &);

        // Operators
        // "*": Operation of SingleParticleOperator on State
        State operator*(const State &) const;

        SingleParticleStatePair operator*(const SingleParticleStatePair &) const;

        ComplexHighRes operatorValue(const State &, const State &) const override;

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

        ComplexHighRes operatorValue(const State &, const State &) const override;

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
    explicit HilbertSpace(SystemScale);

    SingleParticleState createSingleParticleState(const SingleParticleScalarFunction &, const Spin & spin = Spin::None, const std::string & label = "") const;

    ScalarField createScalarField() const;

    ScalarField createScalarField(const std::vector<Complex> &) const;

    ScalarField createScalarField(const SingleParticleScalarFunction &) const;

    Operator createOperator(const SingleParticleFunction &, const SingleParticleFunction &) const;

    Operator createOperator(const DoubleParticleScalarFunction &) const;

    // TODO: proper naming of this function.
    ComplexHighRes operatorValue(const State &, const Operator &, const State &) const;

    double expectationValue(const State &, const Operator &) const;

private:
    SystemScale scale;
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