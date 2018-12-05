#pragma once

#include "hilbert_space.hpp"

struct Setting {
    static Setting defaultSetting();
    
    // NOTE: Parameters are in SI unit system.

    // a:               Width of FD state   (m)
    // d:               Inter-dot distance  (m)
    // B:               B-field             (T)

    SystemScale scale;

    double a, d, B;

    double omegaL() const;

    double omega0() const;

    double omega() const;

    double coulombConstant() const;

    double FDConstant() const;

    double magneticLength() const;
};

enum class Ansatz {
    HL,     // Heitler-London
    HM     // Hund-Mulliken
    // HMp     // Hund-Mulliken + (1+,1) & (1,1+)
};

enum class Basis {
    FD      // Fock-Darwin
};

enum class Orientation {
    Left,
    Right
};

Matrix
coulombMatrix(const Setting &, Ansatz, Basis);

double
calculateJWithSetting_HL(const Setting &);

SingleParticleScalarFunction
fockDarwin(const Setting &, Orientation);

SingleParticleFunction
kineticEnergy(const Setting &setting);

SingleParticleFunction
potentialEnergy(const Setting &setting);

DoubleParticleScalarFunction
coulombEnergy(const Setting &setting);

// Util

SPState
fockDarwinWithSpin(const HilbertSpace &, Setting, Orientation, Spin);

SPState
OrthofockDarwinWithSpin(const HilbertSpace &, Setting, Orientation, Spin);

// For test
DoubleParticleScalarFunction
identity_twoSite(const Setting &setting);