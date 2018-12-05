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

void plotter::plotField(const ScalarField &field, const int trunc, std::string path) {

    int width = field.getScale().width;
    int height = field.getScale().height;
    std::vector<Complex> data = field.getDatas();

    std::ofstream outputFile(path);
    outputFile << std::scientific;

    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j)
            if (((i + j) % trunc == 0) || ((i - j) % trunc == 0))
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