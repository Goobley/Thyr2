#include "Simulation.hpp"
#include <tuple>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "gauss_legendre.c"
#include <cstdio>
#include <xmmintrin.h>
#include <future>

// To enable floating point exception trapping
// _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
// // to disable - will enable if disabled
// _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() ^ _MM_MASK_INVALID);
namespace Constant
{
    constexpr f64 Pi = M_PI;
    constexpr f64 DegToRad = Pi / 180.0;
    /// Astronomical unit [cm]
    constexpr f64 AstroUnit = 1.49597870e13;
    /// arcsec to cm in Sun
    constexpr f64 ArcToCMSun = DegToRad / 3600.0 * AstroUnit;
    /// electron mass [g]
    constexpr f64 ElectronMass = 9.1094e-28;
    /// speed of light [cm/s]
    constexpr f64 SpeedLight = 2.998e10;
    /// electron charge [esu]
    constexpr f64 ElectronCharge = 4.803e-10;
    /// electron rest energy [keV]
    constexpr f64 ElectronRestEnergy = 510.99891;
    /// cgs system to sfu
    constexpr f64 CGSToSFU = 1.0e19;
    /// Boltmann Constant
    constexpr f64 Boltzmann = 1.38e-16;
    /// Solar radius [arcsec]
    constexpr f64 SolarRadius = 960.0;
    /// Planck Constant [erg.s]
    constexpr f64 Planck = 6.6260755e-27;


    /// Gridpoints used in Gauss-Legendre integration
    constexpr int NPTS = 20;
    /// Epsilon convergence precision for floating point in the GyroSim
    constexpr f64 Epsilon = 1e-4;
    /// Gyrofrequency ratio limit to use bessel approximation
    constexpr int BesselLimit = 12;
    /// Max harmonic number before enforcing approx 
    constexpr int BesselHarmonicLimit = 24;
    /// Grid points to evaluate g(Phi)
    constexpr int ND = 100;
}
namespace C = Constant;

template <typename T>
constexpr T
Square(T val)
{
    return val * val;
}

template <typename T>
constexpr T
Cube(T val)
{
    return val * val * val;
}

static f64 Pow(f64 x, f64 y) 
{
    // return std::exp(y * std::log(x));
    return std::pow(x, y);
    // return (T)powf((float)x, (float)y);
}
// #if 0
//     if (y < 0.0)
//       return 1.0 / Pow(x, -y);

//     if (y >= 1.0)
//     {
//         return Square(Pow(x, y / 2.0));
//     } else {
//         f64 low = 0.0;
//         f64 high = 1.0;
//         f64 sqr = std::sqrt(x);
//         f64 acc = sqr;
//         f64 mid = high / 2.0;

//         while (std::abs(mid - y) > tol)
//         {
//              sqr = std::sqrt(sqr);
//              if (mid <= y)
//              {
//                  low = mid;
//                  acc *= sqr;
//              } else {
//                  high = mid;
//                  acc *= (1.0 / sqr);
//              }
//              mid = 0.5 * (low + high);
//         }
//         return acc;
//     }
// #else
//     if (y < 0) return Pow(x, -y);
//     // f64 target = std::pow(x,y);
//     int e = (int)y;
//     union {
//         f64 d;
//         int i[2];
//     } u = { x };
//     u.i[1] = (int)((y - e) * (u.i[1] - 1072632447) + 1072632447);
//     u.i[0] = 0;

//     double r = 1.0;
//     while (e)
//     {
//         if (e & 1)
//         {
//             r *= x;
//         }
//         x *= x;
//         e >>= 1;
//     }
//     f64 val = r * u.d;
//     // f64 err = std::abs(target - val) / target;
//     // if (err > 1e-3)
//         // printf("err %e, %e, %e, x: %e, y: %e\n", err, target, val, x, y);
//     return val;
// #endif

namespace Bessel
{
    f64 Bessj0(f64 x)
    {
        f64 p[5] = {1.e0,-0.1098628627e-2,0.2734510407e-4,
                    -0.2073370639e-5,0.2093887211e-6};
        f64 q[5] = {-0.1562499995e-1,
                    0.1430488765e-3,-0.6911147651e-5,0.7621095161e-6,-0.934945152e-7};
        f64 r[6] = {57568490574.0e0,-13362590354.e0,651619640.7e0,
                    -11214424.18e0,77392.33017e0,-184.9052456e0};
        f64 s[6] = {57568490411.e0,1029532985.0e0,
                    9494680.718e0,59272.64853e0,267.8532712e0,1.0e0};

        if (std::abs(x) < 8.0)
        {
            f64 y = Square(x);
            return (r[0]+y*(r[1]+y*(r[2]+y*(r[3]+y*(r[4]+y*r[5]))))) /
                (s[0]+y*(s[1]+y*(s[2]+y*(s[3]+y*(s[4]+y*s[5])))));
        } else {
            f64 ax = std::abs(x);
            f64 z = 8.0 / ax;
            f64 y = Square(z);
            f64 xx = ax - 0.785398164;

            return sqrt(0.636619772/ax)*(cos(xx)*(p[0]+y*(p[1]+y*(p[2]+y*(p[3]+y*p[4])))) -
                                         z*sin(xx)*(q[0]+y*(q[1]+y*(q[2]+y*(q[3]+y*q[4])))));
        }
    }

    f64 Bessj1(f64 x)
    {
        double r[6] = {72362614232.0e0,-7895059235.0e0,242396853.1e0,
                       -2972611.439e0,15704.48260e0,-30.16036606e0};
        double s[6] = {144725228442.0e0,2300535178.0e0,
                       18583304.74e0,99447.43394e0,376.9991397e0,1.0e0};
        double p[5] = {1.e0,0.183105e-2,-0.3516396496e-4,0.2457520174e-5,
                       -0.240337019e-6};
        double q[5] = {.04687499995e0,-0.2002690873e-3,
                       0.8449199096e-5,-0.88228987e-6,0.105787412e-6};

        if (std::abs(x) < 8.0)
        {
            f64 y = Square(x);

            return x*(r[0]+y*(r[1]+y*(r[2]+y*(r[3]+y*(r[4]+y*r[5]))))) /
                (s[0]+y*(s[1]+y*(s[2]+y*(s[3]+y*(s[4]+y*s[5])))));
        } else {
            f64 ax = std::abs(x);
            f64 z = 8.0 / ax;
            f64 y = Square(z);
            f64 xx = ax - 2.356194491;
            return sqrt(0.636619772/ax)*(cos(xx)*(p[0]+y*(p[1]+y*(p[2]+y*(p[3]+y*p[4])))) -
                                         z*sin(xx)*(q[0]+y*(q[1]+y*(q[2]+y*(q[3]+y*q[4])))))*(x/ax);//*sign(1.,x);
        }
    }

    f64 Bessj(uint n, f64 x)
    {
        const int iacc = 40;
        const double bigno = 1.0e10;
        const double bigni = 1.0e-10;

        f64 tox = 2.0 / x;

        // Should be double, maybe?
        if (x > (f64)n)
        {
            f64 bjm = Bessj0(x);
            f64 bj = Bessj1(x);

            for (uint j = 1; j < n; ++j)
            {
                f64 bjp = j * tox * bj - bjm;
                bjm = bj;
                bj = bjp;
            }

            return bj;
        } else {
            int m = 2 * ((n + (int)(sqrt( (f64)(iacc * n))) )/ 2);
            f64 bessj = 0.0;
            f64 jsum = 0.0;
            f64 sum = 0.0;
            f64 bjp = 0.0;
            f64 bj = 1.0;
            for (int j = m; j >= 1; --j)
            {
                f64 bjm = j * tox * bj - bjp;
                bjp = bj;
                bj = bjm;
                if (std::abs(bj) > bigno)
                {
                    bj *= bigni;
                    bjp *= bigni;
                    bessj *= bigni;
                    sum *= bigni;
                }

                if (jsum != 0)
                    sum += bj;

                jsum = 1 - jsum;

                if (j == (int)n)
                    bessj = bjp;
            }
            sum *= 2.0;
            sum -= bj;
            bessj /= sum;
            return bessj;
        }
    }

    f64 Bessel(uint n, f64 x)
    {
        // BesselMax = std::max(x, BesselMax);
        // BesselMin = std::min(x, BesselMin);

        switch (n)
        {
        case 0:
            return Bessj0(x);
        case 1:
            return Bessj1(x);
        default: // >= 2
            return Bessj(n, x);
        }
    }

    std::pair<f64,f64> BesselWHApprox(uint n, f64 x)
    {
        constexpr f64 A1 = 0.503297;
        constexpr f64 B1 = 1.193000;

        x = x / n;

        f64 a = Pow(Pow(1.0 - Square(x), 1.5) + A1/n, (1.0 / 6.0));
        f64 b = Pow(Pow(1.0 - Square(x), 1.5) + B1/n, (1.0 / 6.0))
            * (1.0 - 1.0 / (5.0 * Pow((f64)n, 2.0 / 3.0)));

        f64 z = (x * exp(sqrt(1.0 - Square(x)))) / (1.0 + sqrt(1.0 - Square(x)));

        // DoubleVec out(2);
        // out.Push(1.0 / sqrt(2.0 * Constant::Pi * n) * pow(z,n) / a);
        // out.Push(a * b * out[0] / x);
        f64 out0 = 1.0 / sqrt(2.0 * Constant::Pi * n) * std::pow(z, n) / a;
        f64 out1 = a * b * out0 / x;

        auto out = std::make_pair(out0, out1);
        return out;
    }

    struct BesselLookup
    {
        BesselLookup()
        {
            printf("Computing Bessel Table\n");
            for (int i = 0; i <= sLimit; ++i)
            {
                vals.emplace_back(std::vector<f64>());
                vals.back().reserve(gridSize);
                for (int j = 0; j < gridSize; ++j)
                {
                    // if (j % gridSize / 10 == 0)
                    //     printf("%f%%\n", 100.0 * (i * gridSize + j)/ ((sLimit+1) * gridSize));
                    // vals.back().emplace_back(Bessel(i, upper / (f64)(gridSize-1.0) * (f64)j));
                    vals.back().emplace_back(jn(i, upper / (f64)(gridSize-1.0) * (f64)j));
                }
            }
            printf("Done\n\n");
        }

        f64 operator()(int s, f64 x) const
        {
            if (s > sLimit)
            {
                printf("s = %d", s);
                assert(false);
            }
            if (x > upper) return Bessel(s, x);

            return vals[s][std::floor(step*x)];
        }

        static constexpr int sLimit = C::BesselHarmonicLimit+1;
        static constexpr int gridSize = 1048576;
        static constexpr f64 upper = 10000.0;
        static constexpr f64 step = (f64)gridSize / upper;
        std::vector<std::vector<f64>> vals;
    };
}

DoubleVec NormaliseEnergyDistribution(DoubleVec& energy, DoubleVec& delta)
{

    // Convert energy in array [keV] to Lorentz factor energy
    for (auto i = energy.begin(); i != energy.end(); ++i)
    {
        *i /= C::ElectronRestEnergy;
        *i += 1.0;
    }

    DoubleVec ce(delta.size());
    for (uint i = 0; i < delta.size(); ++i)
    {
        f64 ce1 = Pow(energy[i] - 1.0  , 1.0 - delta[i]);
        f64 ce2 = Pow(energy[i+1] - 1.0, 1.0 - delta[i]);
        ce[i] = (ce1 - ce2) / (delta[i] - 1.0);
    }

    DoubleVec f(delta.size());
    f[0] = 1.0;
    for (uint i = 1; i < delta.size(); ++i)
        f[i] = Pow(energy[i] - 1.0, delta[i] - delta[i-1]);

    f64 v = 0.0;
    for (uint i = 0; i < f.size(); ++i)
    {
        f[i] *= ce[i];
        v += f[i];
    }
    v = 1.0 / v;

    DoubleVec aNormalised(delta.size());

    aNormalised[0] = v;
    for (uint i = 1; i < delta.size(); ++i)
        aNormalised[i] = aNormalised[i-1] * Pow((energy[i] - 1.0), (delta[i] - delta[i-1]));

    return aNormalised;
}

static void GaussQuadrature(const f64 xMin, const f64 xMax,
                            f64* xCoords, f64* weights,
                            const int numPoints)
// NOTE(Chris): The following comment is verbatim from PJAS source
/*
    Given the lower and upper limits of integration x1 and x2,
    this routine returns arrays x[0..n-1]
    and w[0..n-1] of length n, containing the abscissas and
    weights of the Gauss-Legendre n-point
    quadrature formula.
*/
// TODO(Chris): replace numPoints with vector size?
{
    // Epsilon for convergence
    const f64 Epsilon = 1.0e-14;

    // The roots are symmetric in the interval, so we only need to
    // find half of them
    int midPoint = (numPoints + 1) / 2;
    f64 xMid = 0.5 * (xMax + xMin);
    f64 xHalfLength = 0.5 * (xMax - xMin);

    // Loop over desired roots
    for (int i = 0; i < midPoint; ++i)
    {
        f64 z = std::cos(C::Pi * (i + 0.75)/(numPoints + 0.5));
        f64 prevZ;
        // derivative of the current polynomial in the iteration
        f64 polyDeriv;
        // Starting with this approximation to the ith root, we enter
        // the main refinement loop by Newton's method
        do {
            f64 poly1 = 1.0;
            f64 poly2 = 0.0;

            // Loop over the recurrence relation to get the Legendre
            // polynomial evaluated at z
            for (int j = 0; j < numPoints; ++j)
            {
                f64 prevPoly2 = poly2;
                poly2 = poly1;
                poly1 = ((2.0 * j + 1.0) * z * poly2 - j * prevPoly2) / (j + 1.0);
            }
            // p1 is now the desired Legendre polynomial. We next
            // compute pDeriv, its derivative, by a standard relation
            // also involving p2, the polynomial of one lower order
            polyDeriv = numPoints * (z * poly1 - poly2) / (Square(z) - 1.0);
            prevZ = z;
            z = prevZ - poly1 / polyDeriv; // Looking for convergence by Newton's method
        } while(fabs(z - prevZ) > Epsilon);
        // Scale root to desired interval
        xCoords[i] = xMid - xHalfLength * z;
        // also insert symmetric counterpart
        xCoords[numPoints - i - 1] = xMid + xHalfLength * z;
        // Compute the weight
        weights[i] = 2.0 * xHalfLength / ((1.0 - Square(z)) * Square(polyDeriv));
        // and insert its symmetric counterpart
        weights[numPoints - i - 1] = weights[i];
    }
}

static f64 LerpArray(f64 factor, const DoubleVec& xs, const DoubleVec& ys)
{
    // See if the factor is inside the array (it probably is)
    for (uint i = 0; i < xs.size() - 1; ++i)
    {
        // Find the pair surrounding it and lerp
        if (factor >= xs[i] && factor <= xs[i+1])
        {
            return ys[i] + (factor - xs[i]) * (ys[i+1] - ys[i]) / (xs[i+1] - xs[i]);
        }
    }

    // if it isn't then use an end value
    if (factor >= xs.back())
    {
        return ys.back();
    }
    else
    {
        return ys.front();
    }
}

static DoubleVec LagrangeDeriv(const DoubleVec& xs, const DoubleVec& ys)
{
    // NOTE(Chris): Uses 3 point Lagrange Interpolation to estimate the derivative
    DoubleVec x12(xs.size());
    DoubleVec x01(xs.size());
    DoubleVec x02(xs.size());
    DoubleVec derivs(xs.size());

    for (uint i = 1; i < xs.size()-1; ++i)
    {
        x12[i] = xs[i]   - xs[i+1];
        x01[i] = xs[i-1] - xs[i];
        x02[i] = xs[i-1] - xs[i+1];

        derivs[i] = ys[i-1] * (x12[i] / (x01[i] * x02[i]))
            + ys[i] * (1.0/x12[i] - 1.0/x01[i])
            - ys[i+1] * (x01[i] / (x02[i] * x12[i]));
    }

    // Explicit calculation of first and last points
    derivs[0] = ys[0] * (x01[1] + x02[1]) / (x01[1] * x02[1])
        - ys[1] * x02[1] / (x01[1] * x12[1])
        + ys[2] * x01[1] / (x02[1] * x12[1]);

    uint n = xs.size();
    uint n2 = n - 2;
    derivs[n-1] = -ys[n-3] * x12[n2] / (x01[n2] * x02[n2])
        + ys[n-2] * x02[n2] / (x01[n2] * x12[n2])
        - ys[n-1] * (x02[n2] + x12[n2]) / (x02[n2] * x12[n2]);

    return derivs;
}

std::pair<f64,f64> GyroCore(const CoreConstants& core, const CoreVariables& v, const GyroSimData& d)
{
    static const Bessel::BesselLookup bessel;

    // NOTE(Chris): Core function to calulcate the coefficients
    // Function return values
    f64 G = 0.0;
    f64 H = 0.0;
    // s harmonic limits
    int s1 = std::floor(v.gyroFFB * std::sqrt(1.0 - Square(v.refractionIndex) * Square(core.cosView))) + 1;
    int s2 = std::floor(v.gyroFFB * d.energy.back() * (1.0 + v.refractionIndex * core.betaMax * core.cosView));

    f64 tempDenominator = (1.0 - Square(v.refractionIndex) * Square(core.cosView));

    // h and g functions
    DoubleVec gx(C::NPTS);
    DoubleVec hx(C::NPTS);
    // Convergence flags
    bool hConverged = false;
    bool gConverged = false;

    // Loop of values of s
    for (int s = s1; s < s2; ++s)
    {
        // Energy limits in each s
        f64 lim1 = s/v.gyroFFB;
        f64 lim2 = v.refractionIndex * core.cosView
            * std::sqrt( Square(lim1) + Square(v.refractionIndex) * Square(core.cosView) - 1.0 );

        // Determine lowerG and higerG, limits for the Gauss quadrature
        f64 lowerG = (lim1 - lim2)/tempDenominator;
        if (lowerG < d.energy.front())
            lowerG = d.energy.front();

        f64 higherG = (lim1 + lim2)/tempDenominator;
        if (higherG > d.energy.back())
            higherG = d.energy.back();

        // check order
        if (lowerG >= higherG)
            continue; // increase s and try again

        // DoubleVec gamma;
        // DoubleVec weight;
        f64 gamma[C::NPTS];
        f64 weight[C::NPTS];
        // f64 gamma2[C::NPTS];
        // f64 weight2[C::NPTS];
	    scale_gauss_legendre_table(lowerG, higherG, gamma, weight); 
        // GaussQuadrature(lowerG, higherG, gamma2, weight2, C::NPTS);

        // Loop through the gamma points
        for (int i = 0; i < C::NPTS; ++i)
        {
            // velocity at i
            f64 beta = std::sqrt( Square(gamma[i]) - 1.0)/gamma[i];
            f64 cosPitchAngle = (1.0 - s / v.gyroFFB / gamma[i])/beta/core.cosView/v.refractionIndex;
            f64 sinPitchAngle = std::sqrt(1.0 - Square(cosPitchAngle));
            // https://web.archive.org/web/20161219143458/http://http.developer.nvidia.com/Cg/index_stdlib.html
            // https://stackoverflow.com/questions/3380628/fast-arc-cos-algorithm
            auto Acos = [](f64 x) -> f64 {
                f64 negate = f64(x < 0);
                x = std::abs(x);
                f64 ret = -0.0187293;
                ret = ret * x;
                ret = ret + 0.0742610;
                ret = ret * x;
                ret = ret - 0.2121144;
                ret = ret * x;
                ret = ret + 1.5707288;
                ret = ret * sqrt(1.0-x);
                ret = ret - 2 * negate * ret;
                return negate * 3.14159265358979 + ret;
            };
            // f64 pitchAngle = std::acos(cosPitchAngle); // For pitchAngle in range (0,-pi)
            f64 pitchAngle = Acos(cosPitchAngle); // For pitchAngle in range (0,-pi)

            f64 xs = s * v.refractionIndex * beta * core.sinView * sinPitchAngle /
                (1.0 - v.refractionIndex * beta * core.cosView * cosPitchAngle);

            // Find index of energy array for loop
            // NOTE(Chris): Probably gains to be made here
            uint eIndex = 0;
            for (eIndex = 0; eIndex < d.energy.size()-2; ++eIndex)
            {
                if ((gamma[i] > d.energy[eIndex]) && gamma[i] <= d.energy[eIndex+1])
                    break;
            }

            // Value of pitch angle distribution g(phi) for the current loop
            f64 pitchAngleDistValue = 1.0f;
            // Derivative of previous value wrt phi
            f64 pitchAngleDistValueDeriv = 0.0f;
            switch (d.aniso.type)
            {
            case Anisotropy::None:
            {
                pitchAngleDistValue = 1.0f;
                pitchAngleDistValueDeriv = 0.0f;
            } break;

            case Anisotropy::Gaussian:
            {
                f64 gaussOffset = pitchAngle - d.aniso.Gauss.centre;
                pitchAngleDistValue = std::exp(-0.5 * Square(gaussOffset) / Square(d.aniso.Gauss.width));
                pitchAngleDistValueDeriv = (-gaussOffset) / Square(d.aniso.Gauss.width)
                    * std::exp(-0.5 * Square(gaussOffset) / Square(d.aniso.Gauss.width));
            } break;

            case Anisotropy::Array:
            {
                // Just lerp from precalculated values

                pitchAngleDistValue = LerpArray(pitchAngle, d.aniso.Array.phi, d.aniso.Array.gPhi[eIndex]);
                pitchAngleDistValueDeriv = LerpArray(pitchAngle, d.aniso.Array.phi, d.aniso.Array.dGdPhi[eIndex]);
            } break;

            default:
                assert(false && "Other types of anisotropy have not yet been implemented");

            }

            f64 fe = core.aNormalised[eIndex] * Pow(gamma[i] - 1.0, - d.delta[eIndex]);

            f64 besselDerivative;
            f64 besselValue;
            if (v.gyroFFB < C::BesselLimit && s <= C::BesselHarmonicLimit)
            {
                // f64 besselValueTemp = Bessel::Bessel(s+1, xs);
                besselValue = bessel(s, xs);
                f64 besselValueTemp = bessel(s+1, xs);
                besselDerivative = -besselValueTemp + s / xs * besselValue;
            } else {
                // Wild and Hill (1971) approximation
                auto Jb = Bessel::BesselWHApprox(s, xs);
                besselValue = std::get<0>(Jb);
                besselDerivative = std::get<1>(Jb);
            }

            // Zs function (Klein, 1984), 2nd order term
            f64 zs = -beta * sinPitchAngle * besselDerivative + v.polarisationIndex *
                (core.cosView / core.sinView / v.refractionIndex - beta * cosPitchAngle / core.sinView) * besselValue;

            gx[i] = core.bNormalised * Square(zs) * pitchAngleDistValue *
                (2.0 * C::Pi / core.cosView / (1.0+ Square(v.polarisationIndex))) * v.gyroFFB / beta * fe;

            hx[i] = 0.0;

            if (!hConverged)
            {
                f64 h3a = (d.delta[eIndex] * gamma[i] * (gamma[i] + 1.0) + 2.0 * Square(gamma[i]) - 1.0) /
                    gamma[i] / (Square(gamma[i]) - 1.0);

                f64 tanPitchAngle = (v.refractionIndex * beta * core.cosView - cosPitchAngle) /
                    (gamma[i] * Square(beta) * sinPitchAngle);

                f64 h3b = 0.0;
                if (d.aniso.type != Anisotropy::None)
                    h3b = tanPitchAngle / pitchAngleDistValue * pitchAngleDistValueDeriv;

                // H function
                hx[i] = (h3a + h3b) * gx[i] / Square(v.gyroFFB) / Square(v.refractionIndex);
            }
        } // End iteration over gamma

        // Old value of G
        f64 prevG = G;

        for (int i = 0; i < C::NPTS; ++i)
            G += gx[i] * weight[i];

        f64 gError = std::abs((G - prevG)/prevG);
        if (gError < C::Epsilon)
            gConverged = true;

        if (!hConverged)
        {
            f64 prevH = H;
            for (int i = 0; i < C::NPTS; ++i)
                H += hx[i] * weight[i];

            f64 hError = std::abs((H - prevH) / prevH);
            if (hError < C::Epsilon)
                hConverged = true;
        }

        if (hConverged && gConverged)
            break;


    }
    return std::make_pair(G, H);
}


void CalculatePlasmaIndices(PlasmaIndexData& d, const CoreConstants& c, f64 gyroFFB)
{
    f64 br = Square(c.ratioPlasmaGyro) - Square(gyroFFB);
    br = (Square(c.ratioPlasmaGyro) - Square(gyroFFB));
    f64 ar = std::sqrt(Square(Square(gyroFFB)) * Square(Square(c.sinView)) + 4.0 * Square(gyroFFB) * Square(br) * Square(c.cosView));
    f64 cr = Square(gyroFFB) * Square(c.sinView);

    f64 anum = 2.0 * Square(c.ratioPlasmaGyro) * br;

    f64 dnum1 = +ar - 2.0 * Square(gyroFFB) * br - cr;
    f64 dnum2 = -ar - 2.0 * Square(gyroFFB) * br - cr;

    d.ordinaryRI = std::sqrt(1.0 + anum / dnum1 );
    d.extraordinaryRI = std::sqrt(1.0 + anum / dnum2);

    f64 aknum = 2.0 * gyroFFB * br * c.cosView;
    f64 dknum1 = +ar - cr;
    f64 dknum2 = -ar - cr;
    d.ordinaryPI = -aknum / dknum1;
    d.extraordinaryPI = -aknum / dknum2;
}

ModeParams OrdinaryModeParameters(const CoreConstants& core, const CoreVariables v, 
                                 const GyroSimData& d, const ModeCoeffs coeffs)
{
    ModeParams p;
    if ((v.gyroFFB > core.ratioPlasmaGyro) && (v.refractionIndex > 0.0))
    {
        auto GH = GyroCore(core, v, d);
        p.j = std::get<0>(GH) * coeffs.emissionCoeff;
        p.k = std::get<1>(GH) * coeffs.absorptionCoeff;
    } else {
        p.j = 0.0;
        p.k = 0.0;
    }
    return p;
}

ModeParams ExtraordinaryModeParameters(const CoreConstants& core, const CoreVariables v, 
                                       const GyroSimData& d, const ModeCoeffs coeffs)
{
    ModeParams p;
    if ((v.gyroFFB > (std::sqrt(Square(core.ratioPlasmaGyro) + 0.25) + 0.5)) && (v.refractionIndex > 0.0))
    {
        auto GH = GyroCore(core, v, d);
        p.j = std::get<0>(GH) * coeffs.emissionCoeff;
        p.k = std::get<1>(GH) * coeffs.absorptionCoeff;
    } else {
        p.j = 0.0;
        p.k = 0.0;
    }
    return p;
}

FrequencyAndCoefficients GyroSimulate(GyroSimData& d)
{

    FrequencyAndCoefficients coeffs;
    CoreConstants core;
    switch(d.aniso.type)
    {
    case Anisotropy::None: // Isotropic
    {
        core.bNormalised = 1.0/4.0/C::Pi;
    } break;

    case Anisotropy::Gaussian:
    {
        // DoubleVec tempPhi, weights;
        f64 tempPhi[C::ND];
        f64 weights[C::ND];
        GaussQuadrature(0.0, Constant::Pi, tempPhi, weights, C::ND);

        f64 sum = 0.0;
        for (uint i = 0; i < C::ND; ++i)
        {
            f64 gaussOffset = tempPhi[i] - d.aniso.Gauss.centre;
            f64 currentTGPhi = std::exp(-0.5 * Square(gaussOffset) / Square(d.aniso.Gauss.width));
            sum += currentTGPhi * weights[i];
        }
        core.bNormalised = 1 / (4.0 * sum);
    } break;

    case Anisotropy::Array:
    {
        f64 tempPhi[C::ND];
        f64 weights[C::ND];
        GaussQuadrature(0.0, Constant::Pi, tempPhi, weights, C::ND);

        for (uint e = 0; e < d.energy.size(); ++e)
        {
            f64 sum = 0.0;
            for (uint i = 0; i < C::ND; ++i)
            {
                // TODO(Chris): Change this to support multi energies
                f64 currentTGPhi = LerpArray(tempPhi[i], d.aniso.Array.phi, d.aniso.Array.gPhi[e]);
                sum += currentTGPhi * weights[i];
            }
            core.bNormalised = 1 / (4.0 * sum);

            // Calculate array of derivatives
            // aniso_.Array.dGdPhi.push_back(DoubleVec());
            d.aniso.Array.dGdPhi.emplace_back(LagrangeDeriv(d.aniso.Array.phi, d.aniso.Array.gPhi[e]));


            // Normalise dGdPhi and gPhi
            for (auto&& v : d.aniso.Array.gPhi[e])
                v *= core.bNormalised;

            for (auto&& v : d.aniso.Array.dGdPhi[e])
                v *= core.bNormalised;
        }
        // TODO(Chris): Why is bNormalised_ set to 1 here? -- Because we applied it earlier
        core.bNormalised = 1.0;

    } break;

    default:
        assert(false && "Other anisotropy types are not implemented yet");
    }

    int freqSize = d.frequency.size();
    coeffs.frequency.reserve(freqSize);
    for (const auto& f : d.frequency) 
        coeffs.frequency.emplace_back(f);
    coeffs.ordinaryEmission.reserve(freqSize);
    coeffs.extraordinaryEmission.reserve(freqSize);
    coeffs.ordinaryAbsorption.reserve(freqSize);
    coeffs.extraordinaryAbsorption.reserve(freqSize);

    core.cosView = std::abs(std::cos(d.angle * C::DegToRad));
    core.sinView = std::sqrt(1.0 - Square(core.cosView));

    f64 gyroFreq = 0.5 * C::ElectronCharge * d.bMag / (C::Pi * C::SpeedLight * C::ElectronMass);

    f64 plasmaFreq = C::ElectronCharge * sqrt(d.np /(C::Pi * C::ElectronMass));
    core.ratioPlasmaGyro = plasmaFreq / gyroFreq;

    core.aNormalised = NormaliseEnergyDistribution(d.energy, d.delta);

    // betaMax_ = sqrt(energy_[energyGridSize_ - 1] * energy_[energyGridSize_ - 1] - 1.0) / energy_[energyGridSize_ - 1];
    core.betaMax = std::sqrt(d.energy.back() * d.energy.back() - 1.0) / d.energy.back();

    f64 emissionCoeff = Cube(C::ElectronCharge) / C::ElectronMass /
        Square(C::SpeedLight) * d.bMag * d.nel;

    f64 absorptionCoeff = 4.0 * Square(C::Pi) * C::ElectronCharge / d.bMag * d.nel;

    // Main simulation loop. This is the one to analyse and optimise
    // #pragma omp parallel for default(none) shared(core, coeffs, emissionCoeff, absorptionCoeff, gyroFreq, freqSize, d)

    // Split this out into a function per frequency and then use std;:async to do quick and dirty parallelisation
#if 0
    coeffs.ordinaryEmission.resize(freqSize);
    coeffs.ordinaryAbsorption.resize(freqSize);
    coeffs.extraordinaryEmission.resize(freqSize);
    coeffs.extraordinaryAbsorption.resize(freqSize);
    for (int i=0; i < freqSize; ++i)
    {
        PlasmaIndexData ref;
        CoreVariables c;

        c.gyroFFB = coeffs.frequency[i] / gyroFreq;
        CalculatePlasmaIndices(ref, core, c.gyroFFB);

        // NOTE(Chris): Ordinary mode calulation

        // TODO(Chris): Shouldn't we reserve the length of these
        // vectors first? The growth should probably be constant
        // rather than log(N)... although the allocator time
        // probably doesn't matter in the greater scheme
        if ((c.gyroFFB > core.ratioPlasmaGyro) && (ref.ordinaryRI > 0.0))
        {
            // Pack data
            // CoreVariables core;
            c.refractionIndex = ref.ordinaryRI;
            c.polarisationIndex = ref.ordinaryPI;

            // Compute
            // DoubleVec GHVec = GyroCore(core);

            auto GH = GyroCore(core, c, d);

            // // Store
            // ordinaryEmission_.push_back(GHVec[0] * emissionCoeff_);
            // ordinaryAbsorption_.push_back(GHVec[1] * absorptionCoeff_);
            coeffs.ordinaryEmission[i] = (std::get<0>(GH) * emissionCoeff);
            coeffs.ordinaryAbsorption[i] = (std::get<1>(GH) * absorptionCoeff);
        } else {
            coeffs.ordinaryEmission[i] = (0.0);
            coeffs.ordinaryAbsorption[i] = (0.0);
        }
        // NOTE(Chris): Extraordinary mode calculation
        if ((c.gyroFFB > (std::sqrt(Square(core.ratioPlasmaGyro) + 0.25) + 0.5)) && (ref.extraordinaryRI > 0.0))
        {
            // Calculation method is the same as for ordinary mode
            // CoreVariables core;
            c.refractionIndex = ref.extraordinaryRI;
            c.polarisationIndex = ref.extraordinaryPI;

            // DoubleVec GHVec = GyroCore(core);
            auto GH = GyroCore(core, c, d);
            // extraordinaryEmission_.push_back(GHVec[0] * emissionCoeff_);
            // extraordinaryAbsorption_.push_back(GHVec[1] * absorptionCoeff_);
            coeffs.extraordinaryEmission[i] = (std::get<0>(GH) * emissionCoeff);
            coeffs.extraordinaryAbsorption[i] = (std::get<1>(GH) * absorptionCoeff);
        } else {
            coeffs.extraordinaryEmission[i] = (0.0);
            coeffs.extraordinaryAbsorption[i] = (0.0);
        }
    }
#else
    ModeCoeffs modeCoeffs;
    modeCoeffs.absorptionCoeff = absorptionCoeff;
    modeCoeffs.emissionCoeff = emissionCoeff;
    std::vector<std::future<ModeParams>> futures;
    // Commented sections are for larger jobs, which were... slower...
    // Must be underlying thread pool. Thanks Clang!
    // std::vector<std::future<std::pair<ModeParams, ModeParams>>> futures;
    for (int i = 0; i < freqSize; ++i)
    {
        f64 gyroFFB = coeffs.frequency[i] / gyroFreq;
        PlasmaIndexData ref;
        CalculatePlasmaIndices(ref, core, gyroFFB);

        // futures.emplace_back(std::async(std::launch::async | std::launch::deferred,
        //                                 ModeParameters, core, d, modeCoeffs, gyroFFB, ref));

        CoreVariables ordinaryCore;
        ordinaryCore.gyroFFB = gyroFFB;
        ordinaryCore.polarisationIndex = ref.ordinaryPI;
        ordinaryCore.refractionIndex = ref.ordinaryRI;

        // !!NOTE!! Don't use async|deferred with older gcc libraries. It doesn't launch
        // any other threads. Just have to hit it with all of them unfortunately.
        futures.emplace_back(std::async(std::launch::async | std::launch::deferred,
                                        OrdinaryModeParameters, core, ordinaryCore, d, modeCoeffs));

        CoreVariables extraordinaryCore;
        extraordinaryCore.gyroFFB = gyroFFB;
        extraordinaryCore.polarisationIndex = ref.extraordinaryPI;
        extraordinaryCore.refractionIndex = ref.extraordinaryRI;

        futures.emplace_back(std::async(std::launch::async | std::launch::deferred,
                                        ExtraordinaryModeParameters, core, extraordinaryCore, d, modeCoeffs));

    }

    coeffs.ordinaryEmission.resize(freqSize);
    coeffs.ordinaryAbsorption.resize(freqSize);
    coeffs.extraordinaryEmission.resize(freqSize);
    coeffs.extraordinaryAbsorption.resize(freqSize);

    for (int i = 0; i < freqSize; ++i)
    {
        auto ordinaryParams = futures[2*i].get();
        coeffs.ordinaryEmission[i] = ordinaryParams.j;
        coeffs.ordinaryAbsorption[i] = ordinaryParams.k;

        auto extraParams = futures[2*i + 1].get();
        coeffs.extraordinaryEmission[i] = extraParams.j;
        coeffs.extraordinaryAbsorption[i] = extraParams.k;
        // auto p = futures[i].get();
        // auto oParams = std::get<0>(p);
        // coeffs.ordinaryEmission[i] = oParams.j;
        // coeffs.ordinaryAbsorption[i] = oParams.k;

        // auto xParams = std::get<1>(p);
        // coeffs.extraordinaryEmission[i] = xParams.j;
        // coeffs.extraordinaryAbsorption[i] = xParams.k;
    }
#endif
    coeffs.plasmaDensity = d.np;

    return coeffs;
}


extern "C"
LocalParams GyroSimulateC(GyroSimDataC* data, ThermalRadiationData* td)
{
    GyroSimData d;
    d.bMag = data->bMag;
    d.angle = data->angle;
    d.np = data->np;
    d.nel = data->nel;

    d.frequency.reserve(data->frequencyLen);
    for (int i = 0; i < data->frequencyLen; ++i)
        d.frequency.emplace_back(data->frequency[i]);

    d.energy.reserve(data->energyLen);
    for (int i = 0; i < data->energyLen; ++i)
        d.energy.emplace_back(data->energy[i]);

    d.delta.reserve(data->deltaLen);
    for (int i = 0; i < data->deltaLen; ++i)
        d.delta.emplace_back(data->delta[i]);

    AnisoData aniso;
    switch (data->aniso.type)
    {
        case None: {
            aniso.type = Anisotropy::None;
        } break;
        case Gaussian: {
            aniso.type = Anisotropy::Gaussian;
            aniso.Gauss.centre = data->aniso.Gauss.centre;
            aniso.Gauss.width = data->aniso.Gauss.width;
        } break;
        case Array: {
            aniso.type = Anisotropy::Array;
            int len = data->aniso.Array.phiLen;
            aniso.Array.phi.reserve(len);
            aniso.Array.gPhi.resize(d.delta.size());
            aniso.Array.dGdPhi.resize(d.delta.size());
            for (int i = 0; i < len; ++i)
            {
                aniso.Array.phi.emplace_back(data->aniso.Array.phi[i]);
                aniso.Array.gPhi[i].reserve(len);
                aniso.Array.dGdPhi[i].reserve(len);
            }

            for (int del = 0; del < d.delta.size(); ++del)
            {
                for (int p = 0; p < len; ++p)
                {
                    aniso.Array.gPhi[del].emplace_back(data->aniso.Array.gPhi[del * len + p]);
                    aniso.Array.dGdPhi[del].emplace_back(data->aniso.Array.dGdPhi[del * len + p]);
                }
            }
        } break;
        case Other: {
            aniso.type = Anisotropy::Other;
            assert(false && "Nothing here is implemented...");
        } break;
    }
    d.aniso = aniso;

    auto out = GyroSimulate(d);

    LocalParams result;
    int len = d.frequency.size();
    result.jo = (f64*)calloc(LocalParamsNum * len, sizeof(f64));
    result.ko = result.jo + len;
    result.jx = result.ko + len;
    result.kx = result.jx + len;
    result.jtherm = result.kx + len;
    result.ktherm = result.jtherm + len;

    for (uint i = 0; i < len; ++i)
    {
        f64 t1 = td->temperature < 2.0e5 ? 17.9 + std::log(std::pow(td->temperature, 1.5)) - std::log(d.frequency[i]) :
            24.5 + std::log(td->temperature) - std::log(d.frequency[i]);

        f64 kDulk = 9.78e-3 * d.np * td->protonDensity
            / Square(d.frequency[i])
            / std::pow(td->temperature, 1.5) * t1;
        f64 jDulk = kDulk * Constant::Boltzmann * td->temperature * Square(d.frequency[i])
            / Square(Constant::SpeedLight);

        result.jtherm[i] = 2.0 * jDulk;
        result.ktherm[i] = kDulk;
    }

    // Apply HMinus effects
    for (uint i = 0; i < len; ++i)
    {
        // k for HMinus
        f64 kHm = d.np * td->neutralHDensity / d.frequency[i];
        kHm *= (1.3727e-25 + (4.3748e-10 - 2.5993e-7 / td->temperature) / d.frequency[i]);
        // temperature effect is already "rolled into" Dulk. Apply to k(H-) before summing
        f64 planckFactor = -Constant::Planck * d.frequency[i] / (Constant::Boltzmann * td->temperature);
        kHm *= (1.0 - std::exp(planckFactor));
        result.ktherm[i] += kHm;
    }

    for (int i = 0; i < len; ++i)
    {
        result.ko[i] = out.ordinaryAbsorption[i];
        result.kx[i] = out.extraordinaryAbsorption[i];
        result.jo[i] = out.ordinaryEmission[i];
        result.jx[i] = out.extraordinaryEmission[i];
    }
    return result;
}

typedef std::vector<std::tuple<LocalParams*, std::size_t>> GridInfo;
GridInfo& AllocatedGrids()
{
    static GridInfo vec;
    return vec;
}

extern "C"
LocalParams* AllocateGrid(int xSize, int ySize, int zSize)
{
    auto alloc = (LocalParams*)calloc(xSize * ySize * zSize, sizeof(LocalParams));
    AllocatedGrids().emplace_back(std::make_tuple(alloc, xSize * ySize * zSize));
    return alloc;
}

extern "C"
void FreeGrid(LocalParams* grid)
{
    auto it = std::find_if(AllocatedGrids().begin(), AllocatedGrids().end(), [grid](GridInfo::value_type a)
            {
                return std::get<0>(a) == grid;
            });

    auto ptrCopy = grid;
    if (it != AllocatedGrids().end())
    {
        auto val = *it;
        auto size = std::get<1>(val);
        for (auto i = 0; i < size; ++i)
        {
            free(grid[i].jo);
        }
        free(grid);
    }

    std::remove_if(AllocatedGrids().begin(), AllocatedGrids().end(), [ptrCopy](GridInfo::value_type a)
            {
                return std::get<0>(a) == ptrCopy;
            });
}

extern "C"
int SerializeGrid(const char* filename, LocalParams* grid, int len, int freqLen)
{
    FILE* file = fopen(filename, "wb");
    if (!file)
        return 1;

    for (int i = 0; i < len; ++i)
    {
        if (grid[i].jo)
        {
            fwrite(&i, sizeof(int), 1, file);
            fwrite(grid[i].jo, sizeof(f64), LocalParamsNum * freqLen, file);
        }
    }

    fflush(file);
    fclose(file);
    return 0;
}

extern "C"
LocalParams* DeserializeGrid(const char* filename, int len, int freqLen)
{
    FILE* file = fopen(filename, "rb");
    if (!file)
        return nullptr;

    auto grid = AllocateGrid(len, 1, 1);

    int idx = 0;
    while (fread(&idx, sizeof(int), 1, file) == 1)
    {
        grid[idx].jo = (f64*)calloc(LocalParamsNum * freqLen, sizeof(f64));
        fread(grid[idx].jo, sizeof(f64), LocalParamsNum * freqLen, file);
        grid[idx].ko = grid[idx].jo + freqLen;
        grid[idx].jx = grid[idx].ko + freqLen;
        grid[idx].kx = grid[idx].jx + freqLen;
        grid[idx].jtherm = grid[idx].kx + freqLen;
        grid[idx].ktherm = grid[idx].jtherm + freqLen;
    }

    fclose(file);
    return grid;
}