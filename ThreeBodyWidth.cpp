#include "Pythia8/ParticleData.h"

#include "ThreeBodyWidth.hh"
#include "id.hh"
#include "io.hh"
#include "math.hh"

namespace hnl {

namespace {

const bool debug = false;

}

// The intended use is for numerical integration (a Gaussian quadrature routine is implemented in the base class), root finding, etc, without having to use function pointers (so will probably become obsolete when moving to C++11.)

// Definition of the function to be encapsulated; base class returns 0.

double ThreeBodyWidth::function(std::vector<double> const& integrands) const {
    return function(integrands[0]);
}

bool ThreeBodyWidth::integrate(double& result, double from, double to, double tolerance) const {
    std::vector<double> args(1);
    return integrateGauss(result, 0, from, to, args, tolerance);
}

// Integrate the encapsulated function, f, over argument number iArg, from xLo, to xHi, using Gaussian quadrature, with tolerance tol.
// Return false if precision target could not be reached.
// Adapted from the CERNLIB DGAUSS routine by K.S. Kolbig.

bool ThreeBodyWidth::integrateGauss(double& result, int iArg, double xLo, double xHi, std::vector<double> args, double tol) const {
    result = 0.0;
    if (iArg >= int(args.size())) return false;
    if (xLo >= xHi) return true;

    // 8-point unweighted.
    static double x8[4] = {  0.96028985649753623, 0.79666647741362674, 0.52553240991632899, 0.18343464249564980 };
    static double w8[4] = {  0.10122853629037626, 0.22238103445337447, 0.31370664587788729, 0.36268378337836198 };
    // 16-point unweighted.
    static double x16[8] = { 0.98940093499164993, 0.94457502307323258, 0.86563120238783174, 0.75540440835500303, 0.61787624440264375, 0.45801677765722739, 0.28160355077925891, 0.09501250983763744 };
    static double w16[8] = {  0.027152459411754095, 0.062253523938647893, 0.095158511682492785, 0.12462897125553387, 0.14959598881657673, 0.16915651939500254,  0.18260341504492359, 0.18945061045506850 };

    // Set up integration region.
    double c = 0.001 / std::abs(xHi - xLo);
    double zLo = xLo;
    double zHi = xHi;

    bool nextbin = true;
    while (nextbin) {
        double zMid = 0.5 * (zHi + zLo); // midpoint
        double zDel = 0.5 * (zHi - zLo); // midpoint, relative to zLo
        // Calculate 8-point and 16-point quadratures.
        double s8 = 0.0;
        for (int i = 0; i < 4; i++) {
            double dz = zDel * x8[i];
            args[iArg] = zMid + dz;
            double f1 = function(args);
            args[iArg] = zMid - dz;
            double f2 = function(args);
            s8 += w8[i] * (f1 + f2);
        }
        s8 *= zDel;
        double s16 = 0.0;
        for (int i = 0; i < 8; i++) {
            double dz = zDel * x16[i];
            args[iArg] = zMid + dz;
            double f1 = function(args);
            args[iArg] = zMid - dz;
            double f2 = function(args);
            s16 += w16[i] * (f1 + f2);
        }
        s16 *= zDel;

        // Precision in this bin OK, add to cumulative and go to next.
        if (std::abs(s16 - s8) < tol * (1 + std::abs(s16))) {
            nextbin = true;
            result += s16;
            // Next bin: LO = end of current, HI = end of integration region.
            zLo = zHi;
            zHi = xHi;
            if (zLo == zHi) nextbin = false;

            // Precision in this bin not OK, subdivide.
        } else {
            if (1.0 + c * std::abs(zDel) == 1.0) {
                // Cannot subdivide further at double precision. Fail.
                std::cout << "\n FunctionEncapsulator::integrateGauss(): cannot " << "reach desired tolerance at double precision." << std::endl;
                result = 0.0 ;
                return false;
            }
            zHi = zMid;
            nextbin = true;
        }
    }
    return true;
}

// Solve f(args) = targetValue for argument iArg, on interval from xLo to xHi, using Brent's method, with tolerance tol and a maximum number of iterations maxIter.
// Return false if precision target could not be reached.

bool ThreeBodyWidth::brent(double& solution, double targetValue, int iArg, double xLo, double xHi, std::vector<double> argsIn, double tol, int maxIter) const {

    // Initialize.
    solution = 0.0;

    // Sanity and range checks.
    if (iArg >= int(argsIn.size())) return false;
    if (xLo > xHi) return false;

    std::vector<double> args(argsIn);
    // Evaluate function - targetValue at lower boundary.
    args[iArg] = xLo;
    double f1 = function(args) - targetValue;
    if (std::abs(f1) < tol) {
        solution = xLo;
        return true;
    }
    // Evaluate function - targetValue at upper boundary.
    args[iArg] = xHi;
    double f2 = function(args) - targetValue;
    if (std::abs(f2) < tol) {
        solution = xHi;
        return true;
    }

    // Check if root is bracketed.
    if (f1 * f2 > 0.0) return false;

    // Start searching for root.
    double x1 = xLo;
    double x2 = xHi;
    double x3 = 0.5 * (xLo + xHi);

    int iter = 0;
    while (++iter < maxIter) {
        // Now check at x = x3.
        args[iArg] = x3;
        double f3 = function(args) - targetValue;
        // Check if tolerance on f has been reached.
        if (std::abs(f3) < tol) {
            solution = x3;
            return true;
        }
        // Is root bracketed in lower or upper half?
        if (f1 * f3 < 0.0) xHi = x3;
        else xLo = x3;
        // Check if tolerance on x has been reached.
        if ((xHi - xLo) < tol * (std::abs(xHi) < 1.0 ? xHi : 1.0)) {
            solution = 0.5 * (xLo + xHi);
            return true;
        }

        // Work out next step to take in x.
        double den = (f2 - f1) * (f3 - f1) * (f2 - f3);
        double num = x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 * (f2 - f3)
                     + f1 * x2 * (f3 - f1);
        double dx = xHi - xLo;
        if (den != 0.0) dx = f3 * num / den;

        // First attempt, using gradient
        double x = x3 + dx;

        // If this was too far, just step to the middle
        if ((xHi - x) * (x - xLo) < 0.0) {
            dx = 0.5 * (xHi - xLo);
            x = xLo + dx;
        }
        if (x < x3) {
            x2 = x3;
            f2 = f3;
        } else {
            x1 = x3;
            f1 = f3;
        }
        x3 = x;
    }
    return false; // Maximum number of iterations exceeded.
}

void MesonThreeBodyWidth::set_pointers(Pythia8::ParticleData* particle_data_) {
    if (debug) print("set pointers");
    particle_data = particle_data_;
}

double MesonThreeBodyWidth::get_width(int from_id_, int neutrino_id, int to_id_, int lepton_id_) {
    if (debug) print("get_width", from_id_, "to", to_id_, "with", lepton_id_, "and", neutrino_id);

    id_from = std::abs(from_id_);
    id_to = std::abs(to_id_);
    auto id_lepton = std::abs(lepton_id_);
    auto id_neutrino = std::abs(neutrino_id);

    if (id_from < 100 || id_from > 1000) print("id from", id_from);
    if (id_to < 100 || id_to > 1000) print("id to", id_to);
    if (id_lepton < 10 || id_lepton > 20) print("id lepton", id_lepton);
    if (id_neutrino < 1E5) print("id neutrino", id_neutrino);

    if (id_to > id_from && id_from != 130 && id_to != 211) print("decaying", id_from, "to", id_to);
    if (id_lepton > 20 || id_lepton < 10) print("lepton?", id_lepton);

    mHat = particle_data->m0(id_from);
    m_to = particle_data->m0(id_to);

    auto y_h = m_to / mHat;
    auto y_l = particle_data->m0(id_lepton) / mHat;
    auto y_N = particle_data->m0(id_neutrino) / mHat;

    mr_h = sqr(y_h);
    mr_l = sqr(y_l);
    mr_N = sqr(y_N);

    double width;
    return integrate(width, sqr(y_l + y_N), sqr(1. - y_h), 1e-3 * y_N) ? width : 0.;
}

double lambda(int id_from, int id_to, bool charged) {
    if (is_neutral_Kaon(id_from) && id_to == 211) return charged ? 0.0267 : 0.0117;
    if (id_from == 321 && id_to == 111) return charged ? 0.0277 : 0.0183;
    print("not a valid kaon decay", id_from, id_to, charged);
    return 0;
}

double MesonThreeBodyWidth::K_form_factor(double q2, bool charged) const {
    return 0.970 * (1 + hnl::lambda(id_from, id_to, charged) * q2 / sqr(particle_data->m0(211)));
}

double FF_d_0(int id) {
    if (is_pi(id)) return 0.6117;
    if (is_K(id)) return 0.7647;
    print("not a proper D decay", id);
    return 0;
}

double MesonThreeBodyWidth::z(double q2) const {
    double t_p = sqr(mHat + m_to);
    double t_0 = (mHat + m_to) * sqr(std::sqrt(mHat) - std::sqrt(m_to));
    double first = std::sqrt(t_p - q2);
    double second = std::sqrt(t_p - t_0);
    return (first - second) / (first + second);
}

double c(int id, bool charged) {
    if (is_pi(id)) return charged ? 1.985 : 1.188;
    if (is_K(id)) return charged ? 0.066 : 2.084;
    print("not a proper D decay", id, charged);
    return 0;
}

double P(int id, bool charged) { // GeV^-2
    if (is_pi(id)) return charged ? 0.1314 : 0.0342;
    if (is_K(id)) return charged ? 0.224 : 0;
    print("not a proper D decay", id, charged);
    return 0;
}

double D_s_form_factor(double q2, bool charged) {
    double f_p_eta = 0.495;
    double alpha_p_eta = 0.198;
    double m_D_s = 2.112; // GeV
    double f_0_eta = f_p_eta;
    double alpha_0_eta = 0;
    double q_r = q2 / sqr(m_D_s);
    return charged ? f_p_eta / (1 - q_r) / (1 - alpha_p_eta * q_r) :  f_0_eta / (1 - alpha_0_eta * q_r);
}

double MesonThreeBodyWidth::D_form_factor(double q2, bool charged) const {
    return (FF_d_0(id_to) - c(id_to, charged) * (z(q2) - z(0)) * (1. + (z(q2) + z(0)) / 2.)) / (1. - P(id_to, charged) * q2);
}

double m_pole_prefactor(double q2, int id, bool charged) {
    if (is_D(id) || is_D_s(id)) return 1.;
    if (is_pi(id) || is_K(id)) {
        double xi = q2 / sqr(charged ? 5.325 : 5.65);
        return 1. / (1. - xi);
    }
    print("not a proper B decay", id, charged);
    return 0.;
}

double a_0(int id, bool charged) {
    if (is_pi(id)) return charged ? 0.404 : 0.490;
    if (is_K(id)) return charged ? 0.360 : 0.233;
    if (is_D(id) || is_D_s(id)) return charged ? 0.909 : 0.794;
    print("not a proper B decay", id, charged);
    return 0.;
}

double a_1(int id, bool charged) {
    if (is_pi(id)) return charged ? -0.68 : -1.61;
    if (is_K(id)) return charged ? -0.828 : 0.197;
    if (is_D(id) || is_D_s(id)) return charged ? -7.11 : -2.45;
    print("not a proper B decay", id, charged);
    return 0.;
}

double a_2(int id, bool charged) {
    if (is_pi(id)) return charged ? -0.86 : 0.93;
    if (is_K(id)) return charged ? 1.1 : 0.18;
    if (is_D(id) || is_D_s(id)) return charged ? 66 : 33;
    print("not a proper B decay", id, charged);
    return 0.;
}

double a(int n, int id, int charged) {
    switch (n) {
    case 0 :
        return a_0(id, charged);
    case 1 :
        return a_1(id, charged);
    case 2 :
        return a_2(id, charged);
    default :
        print("Case not covered");
    }
    return 0;
}

double zq2n(int n, double zq2, double zq2N3) {
    switch (n) {
    case 0 :
        return 1.;
    case 1 :
        return zq2 - zq2N3;
    case 2 :
        return sqr(zq2) + 2 * zq2N3;
    default :
        print("unexpected integer");
    }
    return 0.;
}

double MesonThreeBodyWidth::B_form_factor(double q2, bool charged) const {
    double zq2 = z(q2);
    int N = 2;
    double zq2N3 = cube(zq2) / N;
    double sum = 0.;
    for (int n = 0; n < N; ++n) sum += a(n, id_to, charged) * zq2n(n, zq2, zq2N3);
    return m_pole_prefactor(q2, id_to, charged) * sum;
}

double MesonThreeBodyWidth::form_factor(double q2, bool charged) const {
    if (is_K(id_from)) return K_form_factor(q2, charged);
    if (is_D(id_from)) return D_form_factor(q2, charged);
    if (is_D_s(id_from)) return D_s_form_factor(q2, charged);
    if (is_B(id_from) || is_B_s(id_from)) return B_form_factor(q2, charged);
    print("not a proper pseudoscalar decay", id_from, id_to, charged, q2);
    return 0.;
}

double MesonThreeBodyWidth::form_factor_plus(double sqr_q) const {
    return form_factor(sqr_q, true);
}

double MesonThreeBodyWidth::form_factor_0(double sqr_q) const {
    return form_factor(sqr_q, false);
}

double fV(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 1.03;
    if (is_B(id_from) && is_D_star(id_to)) return 0.76;
    if (is_B(id_from) && is_rho(id_to)) return 0.295;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.95;
    if (is_B_s(id_from) && is_K_star(id_to)) return 0.291;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double fA0(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 0.76;
    if (is_B(id_from) && is_D_star(id_to)) return 0.69;
    if (is_B(id_from) && is_rho(id_to)) return 0.231;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.67;
    if (is_B_s(id_from) && is_K_star(id_to)) return 0.289;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double fA1(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 0.66;
    if (is_B(id_from) && is_D_star(id_to)) return 0.66;
    if (is_B(id_from) && is_rho(id_to)) return 0.269;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.70;
    if (is_B_s(id_from) && is_K_star(id_to)) return 0.287;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double fA2(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 0.49;
    if (is_B(id_from) && is_D_star(id_to)) return 0.62;
    if (is_B(id_from) && is_rho(id_to)) return 0.282;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.75;
    if (is_B_s(id_from) && is_K_star(id_to)) return 0.286;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double sigmaV(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 0.27;
    if (is_B(id_from) && is_D_star(id_to)) return 0.57;
    if (is_B(id_from) && is_rho(id_to)) return 0.875;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.372;
    if (is_B_s(id_from) && is_K_star(id_to)) return -0.516;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double sigmaA0(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 0.17;
    if (is_B(id_from) && is_D_star(id_to)) return 0.59;
    if (is_B(id_from) && is_rho(id_to)) return 0.796;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.350;
    if (is_B_s(id_from) && is_K_star(id_to)) return -0.383;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double sigmaA1(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 0.30;
    if (is_B(id_from) && is_D_star(id_to)) return 0.78;
    if (is_B(id_from) && is_rho(id_to)) return 0.54;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.463;
    if (is_B_s(id_from) && is_K_star(id_to)) return 0.;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double sigmaA2(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 0.67;
    if (is_B(id_from) && is_D_star(id_to)) return 1.40;
    if (is_B(id_from) && is_rho(id_to)) return 1.34;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 1.04;
    if (is_B_s(id_from) && is_K_star(id_to)) return 1.05;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double xiV(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 0.;
    if (is_B(id_from) && is_D_star(id_to)) return 0.;
    if (is_B(id_from) && is_rho(id_to)) return 0.;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.561;
    if (is_B_s(id_from) && is_K_star(id_to)) return 2.10;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double xiA0(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 0.;
    if (is_B(id_from) && is_D_star(id_to)) return 0.;
    if (is_B(id_from) && is_rho(id_to)) return 0.055;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.600;
    if (is_B_s(id_from) && is_K_star(id_to)) return 1.58;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double xiA1(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 0.20;
    if (is_B(id_from) && is_D_star(id_to)) return 0.;
    if (is_B(id_from) && is_rho(id_to)) return 0.;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.510;
    if (is_B_s(id_from) && is_K_star(id_to)) return 1.06;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double xiA2(int id_from, int id_to) {
    if (is_D(id_from) && is_K_star(id_to)) return 0.16;
    if (is_B(id_from) && is_D_star(id_to)) return 0.41;
    if (is_B(id_from) && is_rho(id_to)) return -0.21;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.070;
    if (is_B_s(id_from) && is_K_star(id_to)) return -0.074;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double MP(int id_from, int id_to) { // GeV
    if (is_D(id_from) && is_K_star(id_to)) return 1.969;
    if (is_B(id_from) && is_D_star(id_to)) return 6.275;
    if (is_B(id_from) && is_rho(id_to)) return 5.279;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 6.275;
    if (is_B_s(id_from) && is_K_star(id_to)) return 5.367;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double MV(int id_from, int id_to) { // GeV
    if (is_D(id_from) && is_K_star(id_to)) return 2.112;
    if (is_B(id_from) && is_D_star(id_to)) return 6.331;
    if (is_B(id_from) && is_rho(id_to)) return 5.325;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 6.331;
    if (is_B_s(id_from) && is_K_star(id_to)) return 5.415;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double V(double q2, int id_from, int id_to) {
    double ratio = q2 / sqr(MV(id_from, id_to));
    return fV(id_from, id_to) / (1. - ratio) / (1. - sigmaV(id_from, id_to) * ratio - xiV(id_from, id_to) * sqr(ratio));
}

double A0(double q2, int id_from, int id_to) {
    double ratioP = q2 / sqr(MP(id_from, id_to));
    double ratioV = q2 / sqr(MV(id_from, id_to));
    return fA0(id_from, id_to) / (1. - ratioP) / (1. - sigmaA0(id_from, id_to) * ratioV - xiA0(id_from, id_to) * sqr(ratioV));
}

double A1(double q2, int id_from, int id_to) {
    double ratio = q2 / sqr(MV(id_from, id_to));
    return fA1(id_from, id_to) / (1. - sigmaA1(id_from, id_to) * ratio - xiA1(id_from, id_to) * sqr(ratio));
}

double A2(double q2, int id_from, int id_to) {
    double ratio = q2 / sqr(MV(id_from, id_to));
    return fA2(id_from, id_to) / (1. - sigmaA2(id_from, id_to) * ratio - xiA2(id_from, id_to) * sqr(ratio));
}

double MesonThreeBodyWidth::g(double q2) const {
    return V(q2, id_from, id_to) / (mHat + m_to);
}

double MesonThreeBodyWidth::am(double q2) const {
    return (A2(q2, id_from, id_to) * (mHat - m_to) - A1(q2, id_from, id_to) * (mHat - m_to) + 2 * A0(q2, id_from, id_to) * m_to) / q2;
}

double MesonThreeBodyWidth::ff(double q2) const {
    return A1(q2, id_from, id_to) * (mHat + m_to);
}

double MesonThreeBodyWidth::ap(double q2) const {
    return - A2(q2, id_from, id_to) / (mHat + m_to);
}

namespace {

double lambda(double a, double b, double c) {
    return sqr(a) + sqr(b) + sqr(c) - 2 * a * b - 2 * a * c - 2 * b * c;
}

}

double Lambda(double x, double mr_to, double mr_1, double mr_2) {
    return std::sqrt(lambda(1., mr_to, x) * lambda(x, mr_1, mr_2));
}

double MesonThreeBodyWidth::Lambda(double xi) const {
    return hnl::Lambda(xi, mr_h, mr_N, mr_l);
}

double MesonThreeBodyWidth::Gm(double xi) const {
    return xi * (mr_N + mr_l) - sqr(mr_N - mr_l);
}

double MesonThreeBodyWidth::Gp(double xi) const {
    return xi * (mr_N + mr_l) + sqr(mr_N - mr_l);
}

double MesonThreeBodyWidth::F(double xi) const {
    return sqr(1 - xi) - 2 * mr_h * (1 + xi) + sqr(mr_h);
}

double MesonThreeBodyWidth::function(double xi) const {
    double mHat2 = sqr(mHat);
    double q2 = xi * mHat2;
    double Lambdaxi = Lambda(xi);
    double Gmxi = Gm(xi);
    double xi2 = sqr(xi);
    double xi3 = xi2 * xi;
    if (!is_vector(id_to)) {
        double ffp2q2 = sqr(form_factor_plus(q2));
        double term1 = ffp2q2 * cube(Lambdaxi) / 3. / xi3;
        double term2 = ffp2q2 * Lambdaxi * Gmxi * lambda(1., mr_h, xi) / 2. / xi3;
        double term3 = sqr(form_factor_0(q2)) * Lambdaxi * Gmxi * sqr(1. - mr_h) / 2. / xi3;
        return term1 + term2 + term3;
    } else {
        double fq2 = ff(q2);
        double Fxi = F(xi);
        double Gpxi = Gp(xi);
        double apq2 = ap(q2);
        double amq2 = am(q2);
        double term1 = mHat2 * mr_h / 3. / xi2 * g(q2) * Lambdaxi * Fxi * (2 * xi2 - Gpxi);
        double term2 = 1. / 24. / mHat2 / xi3 * sqr(fq2) * Lambdaxi * (3. * Fxi * (xi2 - sqr(mr_l - mr_N)) - sqr(Lambdaxi) + 12. * mr_h * xi * (2 * xi2 - Gpxi));
        double term3 = mHat2 / 24. / xi3 * sqr(apq2) * Lambdaxi * Fxi * (Fxi * (2 * xi2 - Gpxi) + 3. * Gmxi * sqr(1. - mr_h));
        double term4 = mHat2 / 8. / xi * sqr(amq2) * Lambdaxi * Fxi * Gmxi;
        double term5 = 1. / 12. / xi3 * fq2 * apq2 * Lambdaxi * (3 * xi * Fxi * Gmxi + (1 - xi - mr_h) * (3 * Fxi * (xi2 - sqr(mr_l - mr_N)) - sqr(Lambdaxi)));
        double term6 = 1. / 4. / xi2 / fq2 * amq2 * Lambdaxi * Fxi * Gmxi;
        double term7 = mHat2 / 4. / xi2 * apq2 * amq2 * Lambdaxi * Fxi * Gmxi * (1. - mr_h);
        return term1 + term2 + term3 + term4 + term5 + term6 + term7;
    }
}

void NeutrinoThreeBodyWidth::set_pointers(Pythia8::ParticleData* particle_data_, Pythia8::Settings* settings_) {
    if (debug) print("set pointers");
    particle_data = particle_data_;
    settings = settings_;
}

double NeutrinoThreeBodyWidth::get_mass(int id) const {
    switch (id) {
    case 1 :
        return settings->parm("ParticleData:mdRun");
    case 2 :
        return settings->parm("ParticleData:muRun");
    case 3 :
        return settings->parm("ParticleData:msRun");
    case 4 :
        return settings->parm("ParticleData:mcRun");
    case 5 :
        return settings->parm("ParticleData:mbRun");
    case 6 :
        return settings->parm("ParticleData:mtRun");
    default :
        return particle_data->m0(id);
    }
}

double NeutrinoThreeBodyWidth::get_width(int id_from, int id_lepton, int id_up, int id_down) {
    if (debug) print("get_width", id_from, "to", id_up, "with", id_down, "and", id_lepton);

    if (id_from < 1E5) print("id neutrino", id_lepton);
    if (id_lepton < 10 || id_lepton > 20) print("id lepton", id_lepton);
    if (id_down > 20) print("id down", id_down);
    if (id_up > 20) print("id up", id_up);

    auto m_from = get_mass(id_from);

    auto y_l = get_mass(id_lepton) / m_from;
    auto y_u = get_mass(id_up) / m_from;
    auto y_d = get_mass(id_down) / m_from;

    mr1 = sqr(y_l);
    mr2 = sqr(y_d);
    mr3 = sqr(y_u);

    double width;
    return integrate(width, sqr(y_d + y_u), sqr(1. - y_l), 1e-3 * y_l) ? width : 0.;
}

double NeutrinoThreeBodyWidth::function(double x) const {
    return 12. / x * (x - mr2 - mr3) * (1 + mr1 - x) * Lambda(x, mr1, mr2, mr3);
}

}
