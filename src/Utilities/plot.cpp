#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "plot.hpp"
#include "math_utilities.hpp"

void plotter::plotTest() {

}

void plotter::plotFieldAbs(const ScalarField &field) {

}

void plotter::outputToFile(const ScalarField &field, std::string path) {

    int width = field.getWidth();
    int height = field.getHeight();
    std::vector<Complex> data = field.getDatas();

    std::ofstream outputFile(path);
    outputFile << std::scientific;

    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j)
            outputFile << i << " " << j << " " << data[i + j * width].real() << " " << data[i + j * width].imag() << std::endl;

    outputFile.close();
}

void plotter::printMatrix(const Matrix & matrix) {
    for (auto& row : matrix) {
        for (auto& val : row)
            std::cout << std::setprecision(5) << val << " ";
        std::cout << std::endl;
    }
}