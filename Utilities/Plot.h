#pragma once

#include "MathUtilities.h"

namespace plotter {
    void plotTest();
    void plotFieldAbs(const ScalarField &);
    void outputToFile(const ScalarField &, std::string);
}