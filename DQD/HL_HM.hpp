#pragma once

#include "HilbertSpace.hpp"

struct Setting {
    static Setting defaultSetting();
    
    // NOTE: Parameters are in SI unit system.

    // width/height:    System grid.
    // gridSize:        Size of a grid      (m)
    //
    // a:               Width of FD state   (m)
    // d:               Inter-dot distance  (m)
    // B:               B-field             (T)

    int width, height;
    double gridSize;

    double a, d, B;

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

// For test
DoubleParticleScalarFunction
identity_twoSite(const Setting &setting);