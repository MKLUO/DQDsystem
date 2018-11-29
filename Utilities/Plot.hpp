#pragma once

#include "MathUtilities.hpp"

namespace plotter {
    void plotTest();
    void plotFieldAbs(const ScalarField &);
    void outputToFile(const ScalarField &field, std::string path);
}