#include <cmath>

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
ScalarField::operator-(const ScalarField &field) const {
    return *this + field * Complex(-1., 0.);
}

ScalarField
ScalarField::operator*(const SingleParticleScalarFunction &function) const {

    std::vector<Complex> newData;
    newData.resize(data.size());

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            newData[getIndex(x, y)] = function(getX(x), getY(y)) * getData(x, y);

    return ScalarField(width, height, gridSize, newData);
}

ScalarField
ScalarField::operator*(const SingleParticleFunction &function) const {
    return function(*this);
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

ScalarField
ScalarField::operator*(double d) const {
    return *this * Complex(d, 0.);
}

Complex
ScalarField::operator*(const ScalarField &field) const {
    Complex result;

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            result += std::conj(this->getData(x, y)) * field.getData(x, y) * gridSize * gridSize;

    return result;
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

// Some more operators

ScalarField
operator*(Complex c, const ScalarField &field) const {
    return field * c;
}

ScalarField
operator*(double d, const ScalarField &field) const {
    return field * d;
}

ScalarField
operator*(const SingleParticleScalarFunction &function, const ScalarField &field) {
    return field * function;
}

ScalarField
operator*(const SingleParticleFunction &function, const ScalarField &field) const {
    return field * function;
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

// ScalarFields

SingleParticleScalarFunction
        x_field = [](double x, double y) {
    return x;
};

SingleParticleScalarFunction
        y_field = [](double x, double y) {
    return y;
};

SingleParticleScalarFunction
        xx_field = [](double x, double y) {
    return x * x;
};

SingleParticleScalarFunction
        yy_field = [](double x, double y) {
    return y * y;
};

SingleParticleScalarFunction
        sho_field = [](double x, double y) {
    return x * x + y * y;
};

DoubleParticleScalarFunction
        rInv_field = [](double x1, double y1, double x2, double y2) {
    return 1. / hypot((x1 - x2), (y1 - y2));
};

// ScalarFields with settings required

SingleParticleScalarFunction
scalar(Complex c) {
    return [c](double x, double y) {
        return c;
    };
}

SingleParticleScalarFunction
gaussian(double r) {
    return [r](double x, double y) {
        return exp(-(x * x + y * y) / (2. * r * r));
    };
}

SingleParticleScalarFunction
quartic(double a) {
    return [a](double x, double y) {
        return pow((x * x - a * a), 2.0) / (4. * a * a) + y * y;
    };
}

// ScalarFunction

SingleParticleFunction
identity = [](ScalarField field) {
    return field;
};


// Calculate 2D-Laplacian via FFT
// TODO: FFT
SingleParticleFunction
laplacian = [](ScalarField field) {
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
};

// Calculate 2D-"angular momentum" via FFT
SingleParticleFunction
angularMomentum = [](ScalarField field) {
    int width = field.getWidth();
    int height = field.getHeight();
    double gridSize = field.getGridSize();

    std::vector<Complex> field_FT = fourier::fft2d(field.getDatas(), width, height);
    std::vector<Complex> gradX_FT(width * height);
    std::vector<Complex> gradY_FT(width * height);

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y) {
            double kx = double(x) / (2. * M_PI * double(width));
            double ky = double(y) / (2. * M_PI * double(height));

            gradX_FT[x + y * width] = field_FT[x + y * width] * kx * 1.i;
            gradY_FT[x + y * width] = field_FT[x + y * width] * ky * 1.i;
        }

    ScalarField gradX = ScalarField(width, height, gridSize, fourier::ifft2d(gradX_FT, width, height));
    ScalarField gradY = ScalarField(width, height, gridSize, fourier::ifft2d(gradY_FT, width, height));

    return x_field * gradY - y_field * gradX;
};


