#pragma once

#include "HilbertSpace.h"

struct Setting {
    static Setting defaultSetting();

    int width, height;
    double omega0;
};

enum class Orientation {
    Left,
    Right
};

double calculateJWithSetting_HL(const Setting &);

SingleParticleScalarFunction
        fockDarwin(const Setting &, Orientation);

SingleParticleFunction identity();

SingleParticleFunction kineticEnergy(const Setting &setting);

SingleParticleFunction potentialEnergy(const Setting &setting);

DoubleParticleScalarFunction coulombEnergy(const Setting &setting);