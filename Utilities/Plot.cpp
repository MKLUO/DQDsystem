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
void plotter::outputToFile(const ScalarField &field, std::string path) {
    std::ofstream outputFile(path);

    int w = field.getWidth();
    int h = field.getHeight();

    outputFile << std::setprecision(1) << std::scientific;

    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < h; ++j)
            outputFile << field.getData(i, j).real() << "\t";
        outputFile << std::endl;
    }

    outputFile.close();
}