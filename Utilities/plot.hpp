#pragma once

#include "math_utilities.hpp"

namespace plotter {
    void plotTest();
    void plotFieldAbs(const ScalarField &);
    void outputToFile(const ScalarField &field, std::string path);
}