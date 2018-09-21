#pragma once

#include "HilbertSpace.h"

struct Setting
{
    static Setting defaultSetting();

    int width, height;
    double omega0;
};

double calculateJWithSetting_HL(const Setting &);

Complex kineticConstant();

Complex coulombConstant();