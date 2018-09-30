#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "Plot.h"

//#include "gnugraph/GnuGraph.h"

// TODO: implements plotting

void plotter::plotTest() {

}

void plotter::plotFieldAbs(const ScalarField &field) {

}

void plotter::outputToFile(const std::vector<double> &field, int width, int height, std::string path) {
    std::ofstream outputFile(path);

    outputFile << std::scientific;

    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j)
            outputFile << i << " " << j << " " << field[i + j * width] << std::endl;

    outputFile.close();
}