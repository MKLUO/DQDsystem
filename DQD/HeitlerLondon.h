#pragma once

#include "HilbertSpace.h"

struct Setting {
    static Setting defaultSetting();

    // effectiveE:
    // effectiveA:
    // B:
    // alpha: omega0 / omegaL

    // IMPORTANT: Distances are scaled with magneticLength (lb).
    int width, height;
    double gridSize;

    double effectiveE, effectiveA, B, alpha, kappa;

    double omegaL() const;
    double omega0() const;
    double omega() const;
    double coulombConstant() const;
    double FDConstant() const;
    double magneticLength() const;
};

enum class Orientation {
    Left,
    Right
};

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