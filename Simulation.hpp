// Copyright (c) 2015-2018 Christopher Osborne, University of Glasgow
// MIT License: https://opensource.org/licenses/MIT
#include <vector>
typedef double f64;
typedef unsigned int uint;
typedef std::vector<f64> DoubleVec;

struct FrequencyAndCoefficients
{
    f64 plasmaDensity;
    DoubleVec frequency;
    DoubleVec ordinaryAbsorption;
    DoubleVec ordinaryEmission;
    DoubleVec extraordinaryAbsorption;
    DoubleVec extraordinaryEmission;
};


/// Stores type of anisotropy used in the calculations
enum class Anisotropy
{
    None, // Isotropic
    Gaussian, // anisotropic with a gaussian distribution
    Array, // defined by input array distribution
    Other // f(E,U)
};


/// Collection of data related to anisotropy
struct AnisoData
{
    Anisotropy type;

    // NOTE(Chris): I would like to use a union here but the
    // non-PODs prevent it without major headache
    struct
    {
        f64 centre;
        f64 width;
    } Gauss;

    struct
    {
        DoubleVec phi;
        std::vector<DoubleVec> gPhi;
        std::vector<DoubleVec> dGdPhi; // Derivative g wrt phi
    } Array;

    AnisoData(Anisotropy a = Anisotropy::None) : type(a) {};
};

struct CoreConstants
{
    // f64 refractionIndex;
    // f64 polarisationIndex;
    /// Gyro-harmonic
    // f64 gyroFFB;
    f64 betaMax;
    f64 cosView;
    f64 sinView;
    DoubleVec aNormalised;
    f64 bNormalised;
    f64 ratioPlasmaGyro;
};

struct CoreVariables
{
    f64 gyroFFB;
    f64 refractionIndex;
    f64 polarisationIndex;
};

struct PlasmaIndexData
{
    /// Ordinary-mode refractive index
    f64 ordinaryRI;
    /// Extraordinary-moode refractive index
    f64 extraordinaryRI;
    /// Ordinary-mode polarisation index
    f64 ordinaryPI;
    /// Extraordinary-mode polarisation index
    f64 extraordinaryPI;
    /// Gyro-harmonic
    // f64 gyroFFB;
};

struct ModeParams
{
    f64 j;
    f64 k;
};

struct ModeCoeffs
{
    f64 emissionCoeff;
    f64 absorptionCoeff;
};


struct GyroSimData
{
    f64 bMag;
    f64 angle;
    f64 np;
    f64 nel;
    DoubleVec frequency;
    DoubleVec energy;
    DoubleVec delta;
    AnisoData aniso;
};

extern "C"
{
    #include "Simulation.i"
}
// extern "C" typedef struct
// {
//     f64 bMag;
//     f64 angle;
//     f64 np;
//     f64 nel;
//     int frequencyLen;
//     f64* frequency;
//     int energyLen;
//     f64* energy;
//     int deltaLen;
//     f64* delta;
//     // AnisoData aniso;

// } GyroSimDataC;

// extern "C" typedef struct
// {
//     f64 plasmaDensity;
//     int len;
//     f64* frequency;
//     f64* ordinaryAbsorption;
//     f64* ordinaryEmission;
//     f64* extraordinaryAbsorption;
//     f64* extraordinaryEmission;
// } GyroEmissionDataC;
