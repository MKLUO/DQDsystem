#include <vector>
#include <iostream>
#include <string>

using std::vector;
using std::cin;
using std::cout;
using std::string;

#include "Plot.h"

#include "gnugraph/GnuGraph.h"

// TODO: implements plotting

void plotter::plotTest() {
    // TODO: this path should be determined in cmake
    GnuGraph graph("C:/Program Files/gnuplot/bin/gnuplot.exe"); // provide path to executable

    vector<double> x, y;
    for (size_t i = 0; i < 200; ++i)
    {
        x.push_back(double(i));
        y.push_back(sqrt(x[i]));
    }

    const string output = graph.plot(x, y, "y = sqrt(x)");
    cout << output << '\n'; // print any errors to console
    cin.get(); // keep the window open until the user presses ENTER
}
void plotter::plotFieldAbs(const ScalarField &field) {

}