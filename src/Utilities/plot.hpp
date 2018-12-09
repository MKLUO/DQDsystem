#pragma once

#include "math_utilities.hpp"

namespace plotter {
    const std::string tempPath = "temp/TEMPFIELDCAR";

    void plotTest();
    void plotFieldAbs(const ScalarField &);
    void plotField_2D(const ScalarField &field, const int, std::string path);
    void printMatrix(const Matrix &);
}