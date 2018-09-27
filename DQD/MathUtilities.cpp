#include "MathUtilities.h"
#include "Fourier.h"

//////////////////////////////
//      ScalarField        	//
//////////////////////////////

ScalarField::ScalarField(int width_,
                         int height_,
                         double gridSize_,
                         const std::vector<Complex> &data_) {
    width = width_;
    height = height_;
    gridSize = gridSize_;

    data = data_;
}

ScalarField::ScalarField(int width_,
                         int height_,
                         double gridSize_,
                         const SingleParticleScalarFunction &function) {
    width = width_;
    height = height_;
    gridSize = gridSize_;

    data.resize(width * height);

    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j)
            data[getIndex(i, j)] = function(getX(i), getY(j));
}

ScalarField
ScalarField::operator+(const ScalarField &field) const {
    std::vector<Complex> newData;
    newData.resize(data.size());

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            newData[getIndex(x, y)] = this->getData(x, y) + field.getData(x, y);

    return ScalarField(width, height, gridSize, newData);
}

ScalarField
ScalarField::operator*(Complex c) const {
    std::vector<Complex> newData;
    newData.resize(data.size());

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            newData[getIndex(x, y)] = c * this->getData(x, y);

    return ScalarField(width, height, gridSize, newData);
}

Complex
ScalarField::operator*(const ScalarField &field) const {
    Complex result;

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            result += std::conj(this->getData(x, y)) * field.getData(x, y) * gridSize * gridSize;

    return result;
}

ScalarField
ScalarField::operator*(const SingleParticleScalarFunction &function) const {
    std::vector<Complex> newData;
    newData.resize(data.size());

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            newData[getIndex(x, y)] = function(getX(x), getY(y)) * this->getData(x, y);

    return ScalarField(width, height, gridSize, newData);
}

ScalarField
ScalarField::operator*(const SingleParticleFunction &function) const {
    return function(*this);
}

std::vector<Complex>
ScalarField::getDatas() const {
    return data;
}

Complex
ScalarField::getData(int x, int y) const {
    return data[x + y * width];
}

double
ScalarField::getX(int i) const {
    return (i - double(width) / 2) * gridSize;

}

double
ScalarField::getY(int j) const {
    return (j - double(height) / 2) * gridSize;
}

int
ScalarField::getIndex(int i, int j) const {
    return (i + j * width);
}

int
ScalarField::getWidth() const {
    return width;
}

int
ScalarField::getHeight() const {
    return height;
}

double
ScalarField::getGridSize() const {
    return gridSize;
}

Complex &
ScalarField::at(int x, int y) {
    return data[x + y * width];
}


//////////////////////////////
//        Utilities        	//
//////////////////////////////

// TODO: It is supposed to be the most time-consuming part. Consider implement with GPU(CUDA).

Complex
twoSiteIntegral(const ScalarField &left1, const ScalarField &left2,
                const DoubleParticleScalarFunction &function,
                const ScalarField &right1, const ScalarField &right2) {
    Complex result;

    int width = left1.getWidth();
    int height = left1.getHeight();
    double gridSize = left1.getGridSize();

    for (int x1 = 0; x1 < width; ++x1)
        for (int y1 = 0; y1 < height; ++y1)
            for (int x2 = 0; x2 < width; ++x2)
                for (int y2 = 0; y2 < height; ++y2)
                    result += std::conj(left1.getData(x1, y1) * left2.getData(x2, y2)) *
                              function(x1, y1, x2, y2) *
                              right1.getData(x1, y1) * right2.getData(x2, y2);

    return result *
           gridSize * gridSize * gridSize * gridSize;
}


ScalarField
laplacian(const ScalarField &field) {

    //TODO: FFT
    int width = field.getWidth();
    int height = field.getHeight();
    double gridSize = field.getGridSize();

    std::vector<Complex> field_FT = fourier::fft2d(field.getDatas(), width, height);
    std::vector<Complex> result_FT(width * height);

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y) {
            double kx = double(x) / (2. * M_PI * double(width));
            double ky = double(y) / (2. * M_PI * double(height));

            result_FT[x + y * width] = -field_FT[x + y * width] * (kx * kx + ky * ky);
        }

    return ScalarField(width, height, gridSize, fourier::ifft2d(result_FT, width, height));
}

ScalarField
angularMomentum(const ScalarField &) {
    // TODO:
}

Complex
gaussian(double x, double y, double r) {
    return Complex(exp(-(x * x + y * y) / (2. * r * r)), 0.);
}

std::vector<Complex>
fft2d(const std::vector<Complex> &, int width, int height) {

}

std::vector<Complex>
ifft2d(const std::vector<Complex> &, int width, int height) {

}