#pragma once

#include "HilbertSpace.h"

struct Setting {
    static Setting defaultSetting();
    
    // NOTE: Parameters are scaled with dimensionless B=1T quantities.
    // (e.g. magnetic length at 1 Tesla, cyclotron freq. at 1 Tesla)

    // gridSize:    Inter-dot distance / lb(1T)
    // E:           E-field (scaled)
    // a:           Inter-dot distance / lb(1T)
    // B:           B-field / (1T)
    // alpha:       omega0 / omegaL(1T)
    // kappa:       relative dielectric constant

    int width, height;
    double gridSize;

    double E, a, B, alpha, kappa;

    double omegaL() const;

    double omegaL(double) const;

    double omega0() const;

    double omega() const;

    double coulombConstant() const;

    double FDConstant() const;

    double magneticLength() const;

    double magneticLength(double) const;
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