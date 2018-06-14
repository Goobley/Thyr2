// Copyright (c) 2015-2018 Christopher Osborne, University of Glasgow
// MIT License: https://opensource.org/licenses/MIT
typedef double f64;
typedef unsigned int uint;

typedef struct { f64* jo, * ko, * jx, * kx, * jtherm, * ktherm; } LocalParams;
static const int LocalParamsNum = 6;

/// Stores type of anisotropy used in the calculations
enum AnisotropyC
{
    None, // Isotropic
    Gaussian, // anisotropic with a gaussian distribution
    Array, // defined by input array distribution
    Other // f(E,U)
};

/// Collection of data related to anisotropy
struct AnisoDataC
{
    enum AnisotropyC type;

    // NOTE(Chris): I would like to use a union here but the
    // non-PODs prevent it without major headache
    struct
    {
        f64 centre;
        f64 width;
    } Gauss;

    struct
    {
        int phiLen;
        f64* phi; // Angles at which g(phi) is defined
        f64* gPhi; // g is the electron distribution function as a function of angle, 
                   // one array of g(phi) per "chunk" of power law (i.e. number of deltas (spectral indices))
        f64* dGdPhi; // Derivative g wrt phi, same ideas as g(phi)
    } Array;
};

typedef struct
{
    f64 bMag;
    f64 angle;
    f64 np;
    f64 nel;
    int frequencyLen;
    f64* frequency;
    // f64 frequency;
    int energyLen;
    f64* energy;
    int deltaLen;
    f64* delta;
    struct AnisoDataC aniso;

} GyroSimDataC;

typedef struct
{
    f64 temperature;
    f64 neutralHDensity;
    f64 protonDensity;
} ThermalRadiationData;

LocalParams GyroSimulateC(GyroSimDataC* data, ThermalRadiationData* td);
LocalParams* AllocateGrid(int xSize, int ySize, int zSize);
void FreeGrid(LocalParams* grid);
int SerializeGrid(const char* filename, LocalParams* grid, int len, int freqLen);
LocalParams* DeserializeGrid(const char* filename, int len, int freqLen);