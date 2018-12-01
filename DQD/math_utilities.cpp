#define _USE_MATH_DEFINES

#include <math.h>
#include <algorithm>

#include "math_utilities.hpp"
#include "fourier.hpp"

//////////////////////////////
//         Complex          //
//////////////////////////////

ComplexContainer::ComplexContainer() {
    data = std::vector<Complex>({0.0 + 0.0i});
}

ComplexContainer::ComplexContainer(Complex input) {
    data = std::vector<Complex>({input});
}

ComplexContainer::ComplexContainer(double real, double imag) {
    Complex comp = Complex(real, imag);
    data = std::vector<Complex>({comp});
}

ComplexContainer::ComplexContainer(double input) {
    Complex comp = Complex(input);
    data = std::vector<Complex>({comp});
}

ComplexContainer::ComplexContainer(std::vector<Complex> data_) {
    data = data_;
}

ComplexContainer
ComplexContainer::operator+(const ComplexContainer & comp) const {
    
    std::vector<Complex> newData = data;
    if (isZero()) 
        newData = comp.data;
    else 
        newData.insert(newData.end(), comp.data.begin(), comp.data.end());

    ComplexContainer newComp(newData);

    // Shrinking after operation
    if (newComp.size() > COMPLEX_MAX_SIZE) 
        newComp.shrink(COMPLEX_SHRINK_RATIO * double(COMPLEX_MAX_SIZE) / double(newComp.size()));

    return newComp;
}

ComplexContainer
ComplexContainer::operator-(const ComplexContainer & comp) const {
    return *this + (-comp);
}

// When two ComplexContainer multiply and exceed max_size, shrink them inb4.
ComplexContainer
ComplexContainer::operator*(const ComplexContainer & comp) const {
    
    ComplexContainer comp1 = *this;
    ComplexContainer comp2 = comp;

    // Shrinking before operation
    long long oldSize = LONGLONG(comp1.size()) * LONGLONG(comp2.size());
    if (oldSize > LONGLONG(COMPLEX_MAX_SIZE)) {
        double ratio = sqrt(double(COMPLEX_MAX_SIZE) / double(oldSize));
        comp1.shrink(ratio);
        comp2.shrink(ratio);
    }

    std::vector<Complex> newData;
    for (Complex val1 : comp1.data)
        for (Complex val2 : comp2.data)
            newData.push_back(val1 * val2);

    return ComplexContainer(newData);
}

ComplexContainer
ComplexContainer::operator/(const ComplexContainer & comp) const {
   return *this * (1.0 / comp.value());
}

ComplexContainer
ComplexContainer::operator-() const {
    ComplexContainer newComp = *this;
    for (Complex& val : newComp.data)
        val = val * -1.0;

    return newComp;
}

void 
ComplexContainer::operator+=(const ComplexContainer & comp) {
    data = (*this + comp).data;
}

void 
ComplexContainer::operator+=(const Complex & comp) {
    data.push_back(comp);
}


double 
ComplexContainer::real() const {
    return value().real();
}

double 
ComplexContainer::imag() const {
    return value().imag();
}

double
ComplexContainer::norm() const {
    return std::norm(value());
}

ComplexContainer 
ComplexContainer::conj() const {
    ComplexContainer newComp = *this;
    for (Complex& val : newComp.data)
        val = std::conj(val);

    return newComp;
}

int
ComplexContainer::size() const {
    return data.size();
}

bool
ComplexContainer::isZero() const {
    return (data.size() == 1) && (data[0] == 0.0 + 0.0i);
}

Complex 
ComplexContainer::value() const {
    ComplexContainer comp = *this;
    comp.shrink(1);

    return comp.data[0];
}

void 
ComplexContainer::shrink(double ratio) { 
    int newSize = int(std::floor(ratio * double(size())));
    if (newSize < 1) newSize = 1;

    shrink(newSize);
}

void 
ComplexContainer::shrink(int newSize) { 

    if (newSize < 1) newSize = 1;

    //============ (TEST) Method 1 : Naive method ============

    // auto newData = data;

    // std::sort(newData.begin(), newData.end(), 
    //     [](auto val1, auto val2) {
    //         return std::norm(val1) > std::norm(val2);
    //     }
    // );

    // while (newData.size() > newSize) {
    //     newData.end()[-2] = newData.end()[-2] + newData.end()[-1];
    //     newData.pop_back();
    // }
    
    //============ Method 2 : Pairwise reduction ============

    std::vector<double> real_part, imag_part;
    for (Complex val : data) {
        real_part.push_back(val.real());
        imag_part.push_back(val.imag());
    }

    while (real_part.size() > newSize) {
        for (auto newData : {&real_part, &imag_part})
            std::sort(newData->begin(), newData->end(), 
                [](auto val1, auto val2) {
                    return std::abs(val1) > std::abs(val2);
                }
            );

        int reductions;
        if (real_part.size() > 2 * newSize)
            reductions = int(std::floor(0.5 * (double(real_part.size()))));
        else 
            reductions = real_part.size() - newSize;

        std::vector<double> new_real_part, new_imag_part;
        for (int i = 0; i < reductions; i++) {
            new_real_part.push_back(real_part[2*i] + real_part[2*i + 1]);
            new_imag_part.push_back(imag_part[2*i] + imag_part[2*i + 1]);
        }
        for (int i = 2 * reductions; i < real_part.size(); i++) {
            new_real_part.push_back(real_part[i]);
            new_imag_part.push_back(imag_part[i]);
        }

        real_part = new_real_part;
        imag_part = new_imag_part;
    }

    std::vector<Complex> newData;
    for (int i = 0; i < real_part.size(); i++)
        newData.push_back(real_part[i] + 1.0i * imag_part[i]);

    data = newData;
}

ComplexContainer
operator*(double d, const ComplexContainer & comp) {
    return comp * d;
}


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

ComplexContainer
ScalarField::operator*(const ScalarField &field) const {
    ComplexContainer result;

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

ScalarField
ScalarField::conj() const {
    std::vector<Complex> newData(width * height);

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            newData[getIndex(x, y)] = std::conj(this->getData(x, y));

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

// It should be the most time-consuming part. Check the result carefully.
// WARNING: According to the design of Complex object, it might return a Complex which has a data of size width * height. Handle it carefully!
ComplexContainer
twoSiteIntegral(const ScalarField &left1, const ScalarField &left2,
                const DoubleParticleScalarFunction &function,
                const ScalarField &right1, const ScalarField &right2) {

    // Construct Convolution Input

    int width = left1.getWidth();
    int height = left1.getHeight();
    double gridSize = left1.getGridSize();

    ScalarField img(width * 2, height * 2, gridSize);
    ScalarField filter = left1.conj() ^ right1;
    ScalarField weight = left2.conj() ^ right2;

    for (int i = 0; i < 2 * width; ++i)
        for (int j = 0; j < 2 * height; ++j) {
            double x = gridSize * (i - width + 1);
            double y = gridSize * (j - height + 1);
            img.setData(i, j) = function(x, y, 0., 0.);
        }

    ScalarField convResult = fourier::convolution(img, filter);

    convResult = reverse(convResult);

    ComplexContainer result;

    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j)
            result += weight.getData(i, j) * convResult.getData(i, j) * pow(gridSize, 4.0);

    return result;
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
planeWave(double kx, double ky) {
    return [kx, ky](double x, double y) {
        return exp(1.i * (kx * x + ky * y));
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
    double gridSize = field.getGridSize();

    ScalarField field_FT = fourier::fft2d(field);

    ScalarField result_FT(field_FT);

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y) {
            double kx, ky;
            if (x > width / 2 - 1)
                kx = 2. * M_PI * double(x - width) / (double(width) * gridSize);
            else
                kx = 2. * M_PI * double(x) / (double(width) * gridSize);
            if (y > height / 2 - 1)
                ky = 2. * M_PI * double(y - height) / (double(height) * gridSize);
            else
                ky = 2. * M_PI * double(y) / (double(height) * gridSize);

            result_FT.setData(x, y) = field_FT.getData(x, y) * (-1.) * (kx * kx + ky * ky);
        }

    return fourier::ifft2d(result_FT);
};

// Calculate 2D-"angular momentum" via FFT
SingleParticleFunction
        angularMomentum = [](ScalarField field) {

    int width = field.getWidth();
    int height = field.getHeight();
    double gridSize = field.getGridSize();

    ScalarField field_FT = fourier::fft2d(field);

    ScalarField gradX_FT(field_FT);
    ScalarField gradY_FT(field_FT);    

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y) {
            double kx, ky;
            if (x > width / 2 - 1)
                kx = 2. * M_PI * double(x - width) / (double(width) * gridSize);
            else
                kx = 2. * M_PI * double(x) / (double(width) * gridSize);
            if (y > height / 2 - 1)
                ky = 2. * M_PI * double(y - height) / (double(height) * gridSize);
            else
                ky = 2. * M_PI * double(y) / (double(height) * gridSize);

            gradX_FT.setData(x, y) = field_FT.getData(x, y) * kx * Complex(1.i);
            gradY_FT.setData(x, y) = field_FT.getData(x, y) * ky * Complex(1.i);
        }

    ScalarField gradX = fourier::ifft2d(gradX_FT);
    ScalarField gradY = fourier::ifft2d(gradY_FT);

    return x_field * gradY - y_field * gradX;
};


