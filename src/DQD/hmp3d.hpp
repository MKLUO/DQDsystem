#pragma once

#include <fstream>

#include "hilbert_space.hpp"

namespace HMP3D {

    struct Setting {

        Setting(
            const SystemScale, 
            const int, 
            const std::string);
        
        double coulombConstant() const;

        const SystemScale scale;
        const int excitedBasis;
        const std::string infoPath;
    };

    struct EigenState {
        EigenState(
            const ScalarField, 
            const double);
        const ScalarField field;
        const double energy;
    };
    
    std::pair<
        std::vector<EigenState>,
        std::vector<EigenState>>
    parseNN_States(std::string);

    ScalarField
    parseNN_WaveFunction(std::istream);

    Matrix
    hamiltonianMatrix(const Setting &);
    
    DoubleParticleScalarFunction
    coulombEnergy3D(const Setting &);
}