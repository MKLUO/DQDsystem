#pragma once

#include "math_utilities.hpp"

namespace plotter {
    void plotTest();
    void plotFieldAbs(const ScalarField &);
    void plotField(const ScalarField &field, const int, std::string path);
    void printMatrix(const Matrix &);
}