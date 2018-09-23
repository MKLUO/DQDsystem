#pragma once

#include "HilbertSpace.h"

struct Setting
{
    static Setting defaultSetting();

    int width, height;
    double omega0;
};

enum class Orientation
{
    Left,
    Right
};

double calculateJWithSetting_HL(const Setting &);

Complex fockDarwin(double, double, const Setting &, Orientation);

ScalarField kineticEnergy(ScalarField field, const Setting & setting);

ScalarField potentialEnergy(ScalarField field, const Setting & setting);

Complex coulombEnergy(double, double, double, double, const Setting & setting);