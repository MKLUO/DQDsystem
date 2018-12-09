#pragma once

#include "hilbert_space.hpp"

namespace HL_HM { 
    struct Setting {

        Setting(SystemScale, double, double, double);

        static Setting defaultSetting_2D();
        
        // NOTE: Parameters are in SI unit system.

        // a:               Width of FD state   (m)
        // d:               Inter-dot distance  (m)
        // B:               B-field             (T)

        const SystemScale scale;

        double a, d, B;

        double omegaL() const;

        double omega0() const;

        double omega() const;

        double coulombConstant() const;

        double FDConstant() const;

        double magneticLength() const;
    };

    enum class Ansatz {
        HL,     // Heitler-London
        HM     // Hund-Mulliken
    };

    enum class Basis {
        FD      // Fock-Darwin
    };

    enum class Orientation {
        Left,
        Right
    };

    Matrix
    coulombMatrix(const Setting &, Ansatz, Basis);

    double
    calculateJWithSetting_HL(const Setting &);

    SingleParticleScalarFunction
    fockDarwin(const Setting &, Orientation);

    SingleParticleFunction
    kineticEnergy(const Setting &);

    SingleParticleFunction
    potentialEnergy(const Setting &);

    DoubleParticleScalarFunction
    coulombEnergy(const Setting &);

    // Util

    SPState
    fockDarwinWithSpin(const HilbertSpace &, Setting, Orientation, Spin::Type);

    SPState
    OrthofockDarwinWithSpin(const HilbertSpace &, Setting, Orientation, Spin::Type);

    // For test
    DoubleParticleScalarFunction
    identity_twoSite(const Setting &setting);
}