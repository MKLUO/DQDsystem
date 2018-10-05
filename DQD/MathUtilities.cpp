#define _USE_MATH_DEFINES

#include <math.h>

#include "MathUtilities.h"
#include "Fourier.h"

//////////////////////////////
//      ScalarField        	//
//////////////////////////////

ScalarField::ScalarField(int width_,
                         int height_,
                         double gridSize_) {
    width = width_;
    height = height_;
    gridSize = gridSize_;

    data = std::vector<Complex>(width_ * height_);
}

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

ScalarField
ScalarField::operator^(const ScalarField &field) const {
    std::vector<Complex> newData(width * height);

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            newData[getIndex(x, y)] = this->getData(x, y) * field.getData(x, y);

    return ScalarField(width, height, gridSize, newData);
}

std::vector<Complex>
ScalarField::getDatas() const {
    return data;
}

Complex
ScalarField::getData(int x, int y) const {
    return data[x + y * width];
}

Complex &
ScalarField::setData(int x, int y) {
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

std::vector<double>
ScalarField::norm() const {
    std::vector<double> result(width * height);

    for (int i = 0; i < width * height; ++i)
        result[i] = std::norm(data[i]);

    return result;
}

// Some more operators

ScalarField
operator*(Complex c, const ScalarField &field) {
    return field * c;
}

ScalarField
operator*(double d, const ScalarField &field) {
    return field * d;
}

ScalarField
operator*(const SingleParticleScalarFunction &function, const ScalarField &field) {
    return field * function;
}

ScalarField
operator*(const SingleParticleFunction &function, const ScalarField &field) {
    return field * function;
}

//////////////////////////////
//        Utilities        	//
//////////////////////////////

// TODO: It should be the most time-consuming part. Check the result carefully.

Complex
twoSiteIntegral(const ScalarField &left1, const ScalarField &left2,
                const DoubleParticleScalarFunction &function,
                const ScalarField &right1, const ScalarField &right2) {

    // Construct Convolution Input

    int width = left1.getWidth();
    int height = left1.getHeight();
    double gridSize = left1.getGridSize();

    ScalarField img(width * 2, height * 2, gridSize);
    ScalarField filter = left1 ^ right1;
    ScalarField weight = left2 ^ right2;

    for (int i = 0; i < 2 * width; ++i)
        for (int j = 0; j < 2 * height; ++j) {
            double x = gridSize * (i - width + 1);
            double y = gridSize * (j - height + 1);
            img.setData(i, j) = function(x, y, 0., 0.);
        }

    ScalarField convResult = fourier::convolution(img, filter);

    convResult = reverse(convResult);

    Complex result;

    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j) {
            int index = i + j * width;
            result += weight.getData(i, j) * convResult.getData(i, j);
        }
    return result * pow(gridSize, 4.0);
}

ScalarField
reverse(const ScalarField& field) {
    ScalarField result(field);

    int width = field.getWidth();
    int height = field.getHeight();

    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j)
            result.setData(i, j) = field.getData(width - 1 - i,
                                                height - 1 - j);

    return result;
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

DoubleParticleScalarFunction
rInv_field(double gridSize) {
    return [gridSize](double x1, double y1, double x2, double y2) {
        // An approximation is applied around the r=0 point, in which the HilbertSpace gridSize is needed.
        if ((x1 == x2) & (y1 == y2))
            return 2. * sqrt(M_PI) / gridSize;
        else
            return 1. / hypot((x1 - x2), (y1 - y2));
    };
}

// ScalarFunction

SingleParticleFunction
        identity = [](ScalarField field) {
    return field;
};


// Calculate 2D-Laplacian via FFT
SingleParticleFunction
        laplacian = [](ScalarField field) {
    int width = field.getWidth();
    int height = field.getHeight();

    ScalarField field_FT = fourier::fft2d(field);
    ScalarField result_FT(field_FT);

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y) {
            double kx = double(x) / (2. * M_PI * double(width));
            double ky = double(y) / (2. * M_PI * double(height));

            result_FT.setData(x, y) = -field_FT.getData(x, y) * (kx * kx + ky * ky);
        }

    return fourier::ifft2d(result_FT);
};

// Calculate 2D-"angular momentum" via FFT
SingleParticleFunction
        angularMomentum = [](ScalarField field) {

    // TODO: gridSize should be redesigned. It shouldn't be used outside of the HilbertSpace.

    int width = field.getWidth();
    int height = field.getHeight();

    ScalarField field_FT = fourier::fft2d(field);
    ScalarField gradX_FT(field_FT);
    ScalarField gradY_FT(field_FT);

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y) {
            double kx = double(x) / (2. * M_PI * double(width));
            double ky = double(y) / (2. * M_PI * double(height));

            gradX_FT.setData(x, y) = field_FT.getData(x, y) * kx * Complex(1.i);
            gradY_FT.setData(x, y) = field_FT.getData(x, y) * ky * Complex(1.i);
        }

    ScalarField gradX = fourier::ifft2d(gradX_FT);
    ScalarField gradY = fourier::ifft2d(gradY_FT);

    return x_field * gradY - y_field * gradX;
};


