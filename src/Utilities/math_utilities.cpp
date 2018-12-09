#include <math.h>
#include <algorithm>

#include "math_utilities.hpp"
#include "fourier.hpp"

#include "plot.hpp"

//////////////////////////////
//           Grid           //
//////////////////////////////

Grid::Grid(int index_, std::vector<double> coord_):
    index(index_),
    coord(coord_) {}

//////////////////////////////
//       SystemScale        //
//////////////////////////////

// SystemScale::SystemScale(std::vector<int> dims_, double gridSize_):
//     dims(dims_),
//     gridSize(gridSize_),
//     grids(gridding(dims_, gridSize_)) {}

SystemScale::SystemScale(std::vector<int> dims_, double gridSize_):
    dims(padDimToThree(dims_)),
    gridSize(gridSize_) {}

SystemScale
SystemScale::defaultScale_2D() {
    return SystemScale({200, 100}, 1.0E-9);
}

int 
SystemScale::size() const {
    int dimSize = 1;
    for (int dim : dims)
        dimSize *= dim;
    return dimSize;
}

std::vector<double>
SystemScale::coord(int i) const {
    int j = i % (dims[0] * dims[1]);
    double z;
    if (dims[2] == 1)
        z = 0.0;
    else
        z = gridSize * (double((i - j) / (dims[0] * dims[1])) - 0.5 * dims[2]);
    int k = j % (dims[0]);
    double y = gridSize * (double((j - k) / (dims[0])) - 0.5 * dims[1]);
    double x = gridSize * (double(k) - 0.5 * dims[0]);

    return {x, y, z};
}

std::vector<int>
padDimToThree(std::vector<int> dims) {
    if (dims.size() == 2)
        dims.push_back(1);

    if (dims.size() != 3)
        throw std::exception();

    return dims;
}

// std::vector<Grid> *
// gridding(std::vector<int> dims, double gridSize) {
//     std::vector<Grid> * grids = new std::vector<Grid>();    
//     if (dims.size() == 2) {
//         grids->reserve(dims[0] * dims[1]);
//         for (int i = 0; i < dims[0]; ++i)
//         for (int j = 0; j < dims[1]; ++j)
//             grids->emplace_back(
//                 i + j * dims[0],
//                 std::vector<double>
//                 {
//                     (double(i) - dims[0] * 0.5) * gridSize, 
//                     (double(j) - dims[1] * 0.5) * gridSize
//                 }
//             );
//     }
//     else if (dims.size() == 3) {
//         grids->reserve(dims[0] * dims[1] * dims[2]);
//         for (int i = 0; i < dims[0]; ++i)
//         for (int j = 0; j < dims[1]; ++j)
//         for (int k = 0; k < dims[2]; ++k)
//             grids->emplace_back(
//                 i + j * dims[0] + k * dims[0] * dims[1],
//                 std::vector<double>
//                 {
//                     (double(i) - dims[0] * 0.5) * gridSize, 
//                     (double(j) - dims[1] * 0.5) * gridSize,
//                     (double(k) - dims[2] * 0.5) * gridSize
//                 }
//             );
//     }
//     else throw std::exception();

//     return grids;
// }

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

// ComplexContainer::operator Complex() const {
//     return value();
// }

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
    return operator+(-comp);
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
   return operator*(1.0 / comp.value());
}

ComplexContainer
ComplexContainer::operator-() const {
    auto newData = data;
    for (Complex& val : newData)
        val = val * -1.0;

    return ComplexContainer(newData);
}

void 
ComplexContainer::operator+=(const ComplexContainer & comp) {
    data = operator+(comp).data;
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
    auto newData = data;
    for (Complex& val : newData)
        val = std::conj(val);

    return ComplexContainer(newData);
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
                    return std::abs(val1) < std::abs(val2);
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

void 
ComplexContainer::reserve(int newSize) {
    data.reserve(newSize);
}

Complex 
sqrt(const ComplexContainer & c) {
    return sqrt(c.value());
}

ComplexContainer
operator*(double d, const ComplexContainer & comp) {
    return comp * d;
}

std::ostream & 
operator<<(std::ostream & os, const ComplexContainer & comp) {
    os << "(" << comp.real() << ", " << comp.imag() << ")";
    return os;
}

//////////////////////////////
//           Spin           //
//////////////////////////////


Spin::Spin(Type type_) :
    type(type_) {}

double
Spin::operator*(const Spin & spin) const {
    if (*this == spin) 
        return 1.0;
    else if ((type == Type::Up) && (spin.type == Type::Down))
        return 0.0;
    else if ((type == Type::Down) && (spin.type == Type::Up))
        return 0.0;
    else
        throw std::exception();
}

bool 
Spin::operator==(const Spin & spin) const {
    return type == spin.type;
}

Spin::Type
Spin::getType() const {
    return type;
}

std::string
Spin::getLabel() const {
    switch (type) {
        case Spin::Type::Up: return "(+)";
        case Spin::Type::Down: return "(-)";
        case Spin::Type::None: return "";
    }
}

//////////////////////////////
//      ScalarField        	//
//////////////////////////////

ScalarField::ScalarField(const SystemScale scale_) :
    scale(scale_),
    data(std::vector<Complex>(scale_.size())) {}

ScalarField::ScalarField(const SystemScale scale_,
                         const std::vector<Complex> &data_) :
    scale(scale_),
    data(data_) {}

ScalarField::ScalarField(const SystemScale scale_,
                         const SingleParticleScalarFunction &function) :
    scale(scale_),
    data(functionToField(function, scale_)) {}

ScalarField
ScalarField::operator+(const ScalarField &field) const {

    std::vector<Complex> newData(scale.size());

    for (int i = 0; i < scale.size(); ++i)
        newData[i] = data[i] + field.data[i];

    return ScalarField(scale, newData);
}

ScalarField
ScalarField::operator-(const ScalarField &field) const {
    return *this + field * (-1.0);
}

ScalarField
ScalarField::operator*(const SingleParticleScalarFunction &function) const {

    std::vector<Complex> newData(scale.size());

    for (int i = 0; i < scale.size(); ++i)
        newData[i] = function(scale.coord(i)) * data[i];

    return ScalarField(scale, newData);
}

ScalarField
ScalarField::operator*(const SingleParticleFunction &function) const {
    return function(*this);
}

ScalarField
ScalarField::operator*(Complex c) const {

    std::vector<Complex> newData(scale.size());

    for (int i = 0; i < scale.size(); ++i)
        newData[i] = c * data[i];

    return ScalarField(scale, newData);
}

ScalarField
ScalarField::operator/(Complex c) const {
    return operator*(1.0 / c);
}

ScalarField
ScalarField::operator*(double d) const {
    return *this * Complex(d, 0.);
}

ComplexHighRes
ScalarField::operator*(const ScalarField &field) const {
    ComplexHighRes result;

    #ifdef HIGH_RES_COMPLEX
    result.reserve(scale.size());
    #endif

    for (int i = 0; i < scale.size(); ++i)
        result += std::conj(data[i]) * field.data[i] * scale.gridSize * scale.gridSize;

    return result;
}

ScalarField
ScalarField::operator^(const ScalarField &field) const {

    std::vector<Complex> newData(scale.size());

    for (int i = 0; i < scale.size(); ++i)
        newData[i] = data[i] * field.data[i];

    return ScalarField(scale, newData);
}

ScalarField
ScalarField::conj() const {

    std::vector<Complex> newData(scale.size());

    for (int i = 0; i < scale.size(); ++i)
        newData[i] = std::conj(data[i]);

    return ScalarField(scale, newData);
}

std::vector<Complex>
ScalarField::getDatas() const {
    return data;
}

Complex
ScalarField::getData(int i) const {
    return data[i];
}

Complex &
ScalarField::Data(int i) {
    return data[i];
}

SystemScale
ScalarField::getScale() const {
    return scale;
}

std::vector<double>
ScalarField::norm() const {

    std::vector<double> result(scale.size());

    for (int i = 0; i < scale.size(); ++i)
        result[i] = std::norm(data[i]);

    return result;
}

ScalarField
ScalarField::normalize() const {
    ScalarField field = *this;
    return field / sqrt(field * field);
}

void 
ScalarField::plotTemp() const {
    plotter::plotField_2D(*this, 5, plotter::tempPath);
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

std::vector<Complex>
functionToField(
	const SingleParticleScalarFunction & func, 
	const SystemScale scale) {

    std::vector<Complex> newData(scale.size());

    for (int i = 0; i < scale.size(); ++i)
        newData[i] = func(scale.coord(i));

    return newData;
}

// It should be the most time-consuming part. Check the result carefully.
// WARNING: According to the design of Complex object, it might return a Complex which has a data of size width * height. Handle it carefully!
ComplexHighRes
twoSiteIntegral(const ScalarField &left1, const ScalarField &left2,
                const DoubleParticleScalarFunction &function,
                const ScalarField &right1, const ScalarField &right2) {

    // Construct Convolution Input
    // TODO: Check scale match!

    SystemScale oriScale = left1.getScale();
    SystemScale imgScale(
        {oriScale.dims[0] * 2, oriScale.dims[1] * 2},
        oriScale.gridSize);

    ScalarField img(imgScale);
    ScalarField filter = left1.conj() ^ right1;
    ScalarField weight = left2.conj() ^ right2;

    for (int i = 0; i < imgScale.size(); ++i) 
        img.Data(i) = 
            function(
                {
                    imgScale.coord(i)[0] + imgScale.gridSize, 
                    imgScale.coord(i)[1] + imgScale.gridSize,
                    imgScale.coord(i)[2] + imgScale.gridSize
                }, 
                {0., 0., 0.});

    ScalarField convResult = reverse(fourier::convolution2d(img, filter));

    ComplexHighRes result;

    #ifdef HIGH_RES_COMPLEX
    result.reserve(oriScale.size());
    #endif

    for (int i = 0; i < oriScale.size(); ++i) 
        result += 
            weight.getData(i) * 
            convResult.getData(i) * 
            pow(oriScale.gridSize, 4.0);

    return result;
}

ScalarField
reverse(const ScalarField& field) {
    ScalarField result(field);
    int size = field.getScale().size();
    for (int i = 0; i < size; ++i)
        result.Data(i) = 
            field.getData(size - i - 1);

    return result;
}

double 
oneMinus_sqrtOneMinusXX_divideX(const double& x) {
    ComplexHighRes result;
    double lastValue;
    double lastLastValue;
    for (int n = 1; n < 100; n++) {
        double fact1 = 1., fact2 = 1.;
        if (n > 1)
            for (int m = 0; m < n - 1; m++)
                fact1 = fact1 * (m + 2);
        if (n > 2)
            for (int m = 0; m < n - 2; m++)
                fact2 = fact2 * (2.*m + 3);

        double inc = 
            pow(x, 2*n - 1) * fact2 / fact1 / pow(2., n);

        lastLastValue = lastValue;
        lastValue = result.real();
        result += Complex(inc);

        if ((n % 2) && (result.real() == lastLastValue))
            return lastValue;
    }
    return result.real();
}

// Complex
// value(const ComplexHighRes& comp) {
//     return Complex(comp);
// }

// ScalarFields

SingleParticleScalarFunction
        x_field = [](std::vector<double> coord) {
    return coord[0];
};

SingleParticleScalarFunction
        y_field = [](std::vector<double> coord) {
    return coord[1];
};

SingleParticleScalarFunction
        xx_field = [](std::vector<double> coord) {
    return coord[0] * coord[0];
};

SingleParticleScalarFunction
        yy_field = [](std::vector<double> coord) {
    return coord[1] * coord[1];
};

SingleParticleScalarFunction
        sho_field_2D = [](std::vector<double> coord) {
    return xx_field(coord) + yy_field(coord);
};

// ScalarFields with settings required

SingleParticleScalarFunction
scalar(Complex c) {
    return [c](std::vector<double> coord) {
        return c;
    };
}

SingleParticleScalarFunction
planeWave_2D(double kx, double ky) {
    return [kx, ky](std::vector<double> coord) {
        return exp(1.i * (kx * coord[0] + ky * coord[1]));
    };
}

SingleParticleScalarFunction
gaussian_2D(double r) {
    return [r](std::vector<double> coord) {
        return exp(-(sho_field_2D(coord)) / (2. * r * r));
    };
}

SingleParticleScalarFunction
quartic_2D(double a) {
    return [a](std::vector<double> coord) {
        return pow((xx_field(coord) - a * a), 2.0) / (4. * a * a) + yy_field(coord);
    };
}

DoubleParticleScalarFunction
rInv_field_2D(double gridSize) {
    return [gridSize](std::vector<double> coord1, std::vector<double> coord2) {
        // An approximation is applied around the r=0 point, in which the HilbertSpace gridSize is needed.
        // TODO: Verify this approx!
        double x1 = coord1[0];
        double y1 = coord1[1];
        double x2 = coord2[0];
        double y2 = coord2[1];

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
        laplacian_2D = [](ScalarField field) {

    SystemScale scale = field.getScale();

    int widthy = scale.dims[0];
    int height = scale.dims[1];
    double gridSize = scale.gridSize;

    ScalarField field_FT = fourier::fft2d(field);

    ScalarField result_FT(field_FT);

    for (int i = 0; i < scale.size(); ++i) {
        double kx, ky;
        double x = scale.coord(i)[0];
        double y = scale.coord(i)[1];

        if (x >= 0)
            kx =  (x / gridSize - 0.5 * double(widthy)) * 2. * M_PI  / (double(widthy) * gridSize);
        else
            kx =  (x / gridSize + 0.5 * double(widthy)) * 2. * M_PI  / (double(widthy) * gridSize);
        if (y >= 0)
            ky =  (y / gridSize - 0.5 * double(height)) * 2. * M_PI  / (double(height) * gridSize);
        else
            ky =  (y / gridSize + 0.5 * double(height)) * 2. * M_PI  / (double(height) * gridSize);

        result_FT.Data(i) = field_FT.getData(i) * (-1.) * (kx * kx + ky * ky);
    }

    // for (int x = 0; x < width; ++x)
    //     for (int y = 0; y < height; ++y) {
    //         double kx, ky;
    //         if (x > width / 2 - 1)
    //             kx = 2. * M_PI * double(x - width) / (double(width) * gridSize);
    //         else
    //             kx = 2. * M_PI * double(x) / (double(width) * gridSize);
    //         if (y > height / 2 - 1)
    //             ky = 2. * M_PI * double(y - height) / (double(height) * gridSize);
    //         else
    //             ky = 2. * M_PI * double(y) / (double(height) * gridSize);

    //         result_FT.Data(x, y) = field_FT.getData(x, y) * (-1.) * (kx * kx + ky * ky);
    //     }

    return fourier::ifft2d(result_FT);
};

// Calculate 2D-"angular momentum" via FFT
SingleParticleFunction
        angularMomentum_2D = [](ScalarField field) {

    SystemScale scale = field.getScale();

    int widthy = scale.dims[0];
    int height = scale.dims[1];
    double gridSize = scale.gridSize;

    ScalarField field_FT = fourier::fft2d(field);

    ScalarField gradX_FT(field_FT);
    ScalarField gradY_FT(field_FT);

    for (int i = 0; i < scale.size(); ++i) {
        double kx, ky;
        double x = scale.coord(i)[0];
        double y = scale.coord(i)[1];

        if (x >= 0)
            kx =  (x / gridSize - 0.5 * double(widthy)) * 2. * M_PI  / (double(widthy) * gridSize);
        else
            kx =  (x / gridSize + 0.5 * double(widthy)) * 2. * M_PI  / (double(widthy) * gridSize);
        if (y >= 0)
            ky =  (y / gridSize - 0.5 * double(height)) * 2. * M_PI  / (double(height) * gridSize);
        else
            ky =  (y / gridSize + 0.5 * double(height)) * 2. * M_PI  / (double(height) * gridSize);

        gradX_FT.Data(i) = field_FT.getData(i) * kx * Complex(1.i);
        gradY_FT.Data(i) = field_FT.getData(i) * ky * Complex(1.i);
    }

    // for (int x = 0; x < width; ++x)
    //     for (int y = 0; y < height; ++y) {
    //         double kx, ky;
    //         if (x > width / 2 - 1)
    //             kx = 2. * M_PI * double(x - width) / (double(width) * gridSize);
    //         else
    //             kx = 2. * M_PI * double(x) / (double(width) * gridSize);
    //         if (y > height / 2 - 1)
    //             ky = 2. * M_PI * double(y - height) / (double(height) * gridSize);
    //         else
    //             ky = 2. * M_PI * double(y) / (double(height) * gridSize);

    //         gradX_FT.Data(x, y) = field_FT.getData(x, y) * kx * Complex(1.i);
    //         gradY_FT.Data(x, y) = field_FT.getData(x, y) * ky * Complex(1.i);
    //     }

    ScalarField gradX = fourier::ifft2d(gradX_FT);
    ScalarField gradY = fourier::ifft2d(gradY_FT);

    return x_field * gradY - y_field * gradX;
};


