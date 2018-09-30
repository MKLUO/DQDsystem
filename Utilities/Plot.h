#pragma once

#include "MathUtilities.h"

namespace plotter {
    void plotTest();
    void plotFieldAbs(const ScalarField &);
    void outputToFile(const std::vector<double> &, int, int, std::string);
}