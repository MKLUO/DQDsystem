#include "MathUtilities.h"

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

    data.resize(width_ * height_);

    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j) {
            this->at(i, j) = function(getX(i), getY(j));
        }

}

ScalarField
ScalarField::operator+(const ScalarField &field) const {
    std::vector<Complex> newData;
    newData.resize(data.size());

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            newData[getIndex(x, y)] = this->at(x, y) + field.at(x, y);

    return ScalarField(width, height, gridSize, newData);
}

ScalarField
ScalarField::operator*(Complex c) const {
    std::vector<Complex> newData;
    newData.resize(data.size());

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            newData[getIndex(x, y)] = c * this->at(x, y);

    return ScalarField(width, height, gridSize, newData);
}

Complex
ScalarField::operator*(const ScalarField &field) const {
    Complex result;

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            result += std::conj(this->at(x, y)) * field.at(x, y) * gridSize * gridSize;

    return result;
}

ScalarField
ScalarField::operator*(const SingleParticleScalarFunction &function) const {
    std::vector<Complex> newData;
    newData.resize(data.size());

    for (int x = 0; x < width; ++x)
        for (int y = 0; y < height; ++y)
            newData[getIndex(x, y)] = function(getX(x), getY(y)) * this->at(x, y);

    return ScalarField(width, height, gridSize, newData);
}

ScalarField
ScalarField::operator*(const SingleParticleFunction &function) const {
    return function(*this);
}

Complex
ScalarField::at(int x, int y) const {
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

ScalarField
ScalarField::laplacian() const {

}

ScalarField
ScalarField::angularMomentum() const {

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
                    result +=   std::conj(left1.at(x1, y1) * left2.at(x2, y2)) *
                                function(x1, y1, x2, y2) *
                                right1.at(x1, y1) * right2.at(x2, y2) *
                                gridSize * gridSize * gridSize * gridSize;

    return result;
}