#include <cmath>
#include <vector>

#include "FormFactors.hh"
#include "id.hh"
#include "io.hh"
#include "math.hh"

namespace hnl {

namespace {

// Each block below is one fit, with its own masses: as FLAG 2024 (appendix B.3) stresses, the coefficients of a fit
// are only meaningful with the masses, poles and t_0 that went into it. Masses in GeV.

// conformal variable z(q^2, t_0) = (sqrt(t_+ - q^2) - sqrt(t_+ - t_0)) / (sqrt(t_+ - q^2) + sqrt(t_+ - t_0))
double z_map(double q2, double t_plus, double t_0) {
    auto a = std::sqrt(t_plus - q2);
    auto b = std::sqrt(t_plus - t_0);
    return (a - b) / (a + b);
}

// optimal t_0 = t_+ (1 - sqrt(1 - t_- / t_+))
double t_0_optimal(double t_plus, double t_minus) {
    return t_plus * (1. - std::sqrt(1. - t_minus / t_plus));
}

double pole(double q2, double mass) {
    return mass > 0. ? 1. / (1. - q2 / sqr(mass)) : 1.;
}

// BCL parametrisation, eqs. (535)-(545) of FLAG 2024 (arXiv:2411.04268):
// f_+ = 1 / (1 - q^2/m_+^2) sum_n a+_n [z^n - (-1)^(n - N) n / N z^N],  f_0 = 1 / (1 - q^2/m_0^2) sum_n a0_n z^n.
// A pole mass of 0 means no pole. If a0 is one coefficient short, the last one follows from f_+(0) = f_0(0).
struct BCL {
    BCL(double t_plus_, double t_0_, double m_plus_, double m_0_, std::vector<double> a_plus_, std::vector<double> a_0_) :
        t_plus(t_plus_), t_0(t_0_), m_plus(m_plus_), m_0(m_0_), a_plus(std::move(a_plus_)), a_0(std::move(a_0_)) {
        if (a_0.size() == a_plus.size()) return;
        auto z0 = z(0.);
        auto sum = 0.;
        for (auto n = 0u; n < a_0.size(); ++n) sum += a_0[n] * std::pow(z0, n);
        a_0.push_back((f_plus(0.) - sum) / std::pow(z0, a_0.size()));
    }
    double z(double q2) const {
        return z_map(q2, t_plus, t_0);
    }
    double f_plus(double q2) const {
        auto N = static_cast<int>(a_plus.size());
        auto zq = z(q2);
        auto sum = 0.;
        for (auto n = 0; n < N; ++n) sum += a_plus[n] * (std::pow(zq, n) - ((n - N) % 2 == 0 ? 1. : -1.) * n / N * std::pow(zq, N));
        return pole(q2, m_plus) * sum;
    }
    double f_0(double q2) const {
        auto zq = z(q2);
        auto sum = 0.;
        for (auto n = 0u; n < a_0.size(); ++n) sum += a_0[n] * std::pow(zq, n);
        return pole(q2, m_0) * sum;
    }
    double t_plus;
    double t_0;
    double m_plus;
    double m_0;
    std::vector<double> a_plus;
    std::vector<double> a_0;
};

// --- K -> pi: f_+(0) of FLAG 2024 (arXiv:2411.04268, eq. 72, N_f = 2+1+1) for K0 -> pi-. Shape: quadratic f_+ and
//     linear f_0 in t / m_pi+^2 (eq. 28 of FlaviaNet, arXiv:1005.2323) with their average of K_e3 and K_mu3 without
//     NA48 K_mu3 (table 10): lambda'_+ = 25.04e-3, lambda''_+ = 1.57e-3, lambda_0 = 15.90e-3. K+ -> pi0 has the
//     isospin-breaking correction 1 + delta_SU(2), delta_SU(2) = 0.029 (eq. 26 there); the 1/sqrt(2) of the pi0 is
//     in the width (clebsch_gordan_2).
std::array<double, 2> K_pi(int id_from, double q2) {
    auto t = q2 / sqr(0.13957039); // m_pi+, PDG
    auto f_0 = 0.9698 * (id_from == 321 ? 1.029 : 1.);
    return {f_0 * (1. + 25.04e-3 * t + 1.57e-3 / 2. * sqr(t)), f_0 * (1. + 15.90e-3 * t)};
}

// --- D -> K: FLAG 2024 average of ETM 17D, HPQCD 21A and FNAL/MILC 22, table 30 and appendix B.3.1 (N+ = N0 = 4)
BCL const& D_K() {
    static auto const m_D = 1.87265, m_K = 0.495644;
    static auto const t_plus = sqr(m_D + m_K);
    static BCL const bcl(t_plus, t_0_optimal(t_plus, sqr(m_D - m_K)), 2.1122, 2.317, {0.7953, -1.0090, 0.22, 0.14}, {0.7026, 0.773, 0.54});
    return bcl;
}

// --- D -> pi: FNAL/MILC 22 (arXiv:2212.12648), tables X, XII and XIII: t_0 = 0, poles D* and D0*(2300). FLAG 2024
//     quotes no shape for D -> pi because ETM 17D disagrees at high q^2; FNAL/MILC 22 is the more precise calculation.
BCL const& D_pi() {
    static BCL const bcl(sqr(1.86483 + 0.1349768), 0., 2.00685, 2.300, {0.6300, -0.610, -0.20, 0.30}, {0.6300, 0.330, -0.31, -1.90});
    return bcl;
}

// --- D_s -> eta: light-cone sum rule shape of 1805.08567 (table 13: f_+(0) = 0.495, alpha = 0.198, pole m_D_s*,
//     f_0 = f_+(0)), normalised to B(D_s+ -> eta e+ nu) = 2.27(6) % (PDG 2025) with tau_D_s and |V_cs| of PDG 2025.
std::array<double, 2> D_s_eta(double q2) {
    auto const norm = 0.9716; // sqrt(2.27 % / 2.405 %), the latter from this shape
    auto f_0 = 0.495 * norm;
    auto q_r = q2 / sqr(2.112);
    return {f_0 / (1. - q_r) / (1. - 0.198 * q_r), f_0};
}

// --- B -> pi: FLAG 2024 fit to lattice (FNAL/MILC 15, RBC/UKQCD 15, JLQCD 22) and BaBar and Belle data, table 50
//     and appendix B.3.2: f_+ with the B* pole, f_0 without pole
BCL const& B_pi() {
    static auto const m_B = 5.27934, m_pi = 0.1349768;
    static BCL const bcl(sqr(m_B + m_pi), (m_B + m_pi) * sqr(std::sqrt(m_B) - std::sqrt(m_pi)), 5.32471, 0., {0.425, -0.441, -0.52}, {0.560, -1.346});
    return bcl;
}

// --- B_s -> K: FLAG 2024 lattice average, table 40 and appendix B.3.3 (N+ = N0 = 4; the cut starts at B pi)
BCL const& B_s_K() {
    static auto const m_B = 5.27931, m_pi = 0.1349766, m_B_s = 5.36688, m_K = 0.493677;
    static auto const t_plus = sqr(m_B + m_pi);
    static BCL const bcl(t_plus, t_0_optimal(t_plus, sqr(m_B_s - m_K)), 5.3251, 5.68, {0.370, -0.68, 0.55, 2.11}, {0.234, 0.135, 0.20});
    return bcl;
}

// --- B -> D: FLAG 2024 fit to lattice (FNAL/MILC 15C, HPQCD 15) and BaBar and Belle data, table 52 and appendix B.3.5
BCL const& B_D() {
    static auto const m_B = 5.27931, m_D = (1.86483 + 1.86965) / 2.;
    static BCL const bcl(sqr(m_B + m_D), (m_B + m_D) * sqr(std::sqrt(m_B) - std::sqrt(m_D)), 0., 0., {0.8946, -8.03, 50.1}, {0.7804, -3.38});
    return bcl;
}

// --- B_s -> D_s: HPQCD 19 (arXiv:1906.00701), appendix A and table VIII (as reproduced in table 46 of FLAG 2024):
//     t_0 = 0, poles B_c* and B_c0
BCL const& B_s_D_s() {
    static BCL const bcl(sqr(5.3669 + 1.9683), 0., 6.329, 6.704, {0.666, -3.24, -0.1}, {0.666, -0.26, -0.1});
    return bcl;
}

// --- B_c -> eta_c: Ivanov, Korner, Santorelli, hep-ph/0007169, eq. (57) with f(0) of table II and m_fit, delta of
//     table VIII: f(0) / (1 - q^2/m^2 - delta q^4/m^4); f_0 = f_+ + q^2 / (m_P^2 - m_P'^2) f_-. No lattice result
//     is published for the full q^2 range.
double IKS(double q2, double f0, double m, double delta) {
    auto ratio = q2 / sqr(m);
    return f0 / (1. - ratio - delta * sqr(ratio));
}

std::array<double, 2> B_c_eta_c(double q2, double m_from, double m_to) {
    auto f_plus = IKS(q2, 0.76, 6.37, 0.087);
    auto f_minus = IKS(q2, -0.38, 6.22, 0.060);
    return {f_plus, f_plus + q2 / (sqr(m_from) - sqr(m_to)) * f_minus};
}

// --- vector daughters

// Melikhov, Stech, hep-ph/0001113, table 3 (as in table 14 of 1805.08567): f(0), sigma, xi and the poles M_P, M_V
struct MelikhovStech {
    double fV, fA0, fA1, fA2, sigmaV, sigmaA0, sigmaA1, sigmaA2, xiV, xiA0, xiA1, xiA2, MP, MV;
    std::array<double, 4> operator()(double q2) const {
        auto ratioP = q2 / sqr(MP);
        auto ratioV = q2 / sqr(MV);
        return {fV / (1. - ratioV) / (1. - sigmaV * ratioV - xiV * sqr(ratioV)),
                fA0 / (1. - ratioP) / (1. - sigmaA0 * ratioV - xiA0 * sqr(ratioV)),
                fA1 / (1. - sigmaA1 * ratioV - xiA1 * sqr(ratioV)),
                fA2 / (1. - sigmaA2 * ratioV - xiA2 * sqr(ratioV))};
    }
};

// --- D -> K*: Melikhov-Stech shape, normalised to B(D+ -> anti-K*0 e+ nu) = 5.40(10) % (PDG 2025) with tau_D+ and
//     |V_cs| of PDG 2025. There is no lattice calculation; the shape agrees with the measured ratios
//     r_V = V(0)/A_1(0) = 1.49 and r_2 = A_2(0)/A_1(0) = 0.80 (average in table 5 of BESIII, arXiv:2408.04422) to 5 %.
std::array<double, 4> D_K_star(double q2) {
    auto const norm = 0.9417; // sqrt(5.40 % / 6.089 %), the latter from this shape
    static MelikhovStech const ms{1.03, 0.76, 0.66, 0.49, 0.27, 0.17, 0.30, 0.67, 0., 0., 0.20, 0.16, 1.969, 2.112};
    auto ff = ms(q2);
    for (auto& f : ff) f *= norm;
    return ff;
}

// --- B -> rho, B_s -> K*: light-cone sum rules of Bharucha, Straub, Zwicky (BSZ, arXiv:1503.05534), eqs. (15)-(16):
//     F = 1 / (1 - q^2/m_R^2) sum_k alpha_k [z(q^2) - z(0)]^k with t_+- = (m_B +- m_V)^2 and the optimal t_0;
//     poles of table 3 (b -> d column, the resonances of the b u current). They give A_12 instead of A_2.
struct BSZ {
    double m_B, m_V;
    std::array<double, 3> A0, A1, A12, V;
    double F(std::array<double, 3> const& alpha, double q2, double m_R) const {
        auto t_plus = sqr(m_B + m_V);
        auto t_0 = t_0_optimal(t_plus, sqr(m_B - m_V));
        auto dz = z_map(q2, t_plus, t_0) - z_map(0., t_plus, t_0);
        return pole(q2, m_R) * (alpha[0] + alpha[1] * dz + alpha[2] * sqr(dz));
    }
    std::array<double, 4> operator()(double q2) const {
        auto a1 = F(A1, q2, 5.724);
        auto a12 = F(A12, q2, 5.724);
        // A_12 = ((m_B + m_V)^2 (m_B^2 - m_V^2 - q^2) A_1 - lambda A_2) / (16 m_B m_V^2 (m_B + m_V)), eq. (4) of BSZ
        auto lambda = (sqr(m_B + m_V) - q2) * (sqr(m_B - m_V) - q2);
        auto a2 = (sqr(m_B + m_V) * (sqr(m_B) - sqr(m_V) - q2) * a1 - 16. * m_B * sqr(m_V) * (m_B + m_V) * a12) / lambda;
        return {F(V, q2, 5.325), F(A0, q2, 5.279), a1, a2};
    }
};

// B -> rho: BSZ LCSR-only fit, table 14
BSZ const& B_rho() {
    static BSZ const bsz{5.27963, 0.77526, {0.36, -0.83, 1.33}, {0.26, 0.39, 0.16}, {0.30, 0.76, 0.46}, {0.33, -0.86, 1.80}};
    return bsz;
}

// B_s -> K*: BSZ combined LCSR and lattice fit, table 15
BSZ const& B_s_K_star() {
    static BSZ const bsz{5.36689, 0.89167, {0.36, -0.36, 8.03}, {0.22, 0.24, 1.77}, {0.27, 1.12, 3.43}, {0.28, -0.82, 5.08}};
    return bsz;
}

// --- B -> D*: FLAG 2024 BGL fit to lattice (N_f = 2+1: FNAL/MILC 21, JLQCD 23) and Belle and BaBar data, table 53,
//     in the parametrisation of FNAL/MILC 21 (arXiv:2105.14019), eqs. (55)-(70) and tables 9 and 10 (appendix B.3.7
//     of FLAG 2024). g = 2 V / (m_B + m_V), f = (m_B + m_V) A_1, F_1 = 8 m_B m_V A_12 (the helicity combination of A_1
//     and A_2), F_2 = 2 A_0. The coefficients a^F1_0 and a^F2_0 follow from F_1(w = 1) = m_B (1 - r) f(w = 1) and
//     A_0(0) = A_3(0).
struct BGL {
    double m_B = 5.280, m_V = 2.010; // table 10 of FNAL/MILC 21
    double r = m_V / m_B;
    double n_I = 2.6, chi_T_1p = 3.894e-4, chi_T_1m = 5.131e-4, chi_L_1p = 1.9421e-2;
    std::vector<double> poles_1m{6.329, 6.920, 7.020}, poles_1p{6.739, 6.750, 7.145, 7.150}, poles_0m{6.275, 6.842, 7.250};
    std::vector<double> a_g{0.03036, -0.083}, a_f{0.01213, 0.0234, -0.59}, a_F1{0., 0.00141, -0.005}, a_F2{0., -0.093};
    BGL() {
        a_F1[0] = m_B * (1. - r) * a_f[0] * phi_F1(0.) / phi_f(0.);
        // F_2(q^2 = 0) = 2 A_0(0) = 2 A_3(0) = 2 F_1(0) / (m_B^2 - m_V^2)
        auto z_max = z(0.);
        auto F_2 = 2. * F_1(0.) / (sqr(m_B) - sqr(m_V));
        a_F2[0] = F_2 * blaschke(z_max, poles_0m) * phi_F2(z_max) - series(a_F2, z_max);
    }
    double w(double q2) const {
        return (sqr(m_B) + sqr(m_V) - q2) / (2. * m_B * m_V);
    }
    double z_of_w(double w) const {
        return (std::sqrt(w + 1.) - std::sqrt(2.)) / (std::sqrt(w + 1.) + std::sqrt(2.));
    }
    double z(double q2) const {
        return z_of_w(w(q2));
    }
    double blaschke(double z, std::vector<double> const& poles) const {
        auto product = 1.;
        for (auto mass : poles) {
            auto z_p = z_of_w(w(sqr(mass)));
            product *= (z - z_p) / (1. - z * z_p);
        }
        return product;
    }
    double denominator(double z) const {
        return (1. + r) * (1. - z) + 2. * std::sqrt(r) * (1. + z);
    }
    double phi_g(double z) const {
        return 16. * sqr(r) * std::sqrt(n_I / (3. * M_PI * chi_T_1m)) * sqr(1. + z) / std::sqrt(1. - z) / std::pow(denominator(z), 4);
    }
    double phi_f(double z) const {
        return 4. * r / sqr(m_B) * std::sqrt(n_I / (3. * M_PI * chi_T_1p)) * (1. + z) * std::pow(1. - z, 1.5) / std::pow(denominator(z), 4);
    }
    double phi_F1(double z) const {
        return 4. * r / cube(m_B) * std::sqrt(n_I / (6. * M_PI * chi_T_1p)) * (1. + z) * std::pow(1. - z, 2.5) / std::pow(denominator(z), 5);
    }
    double phi_F2(double z) const {
        return 8. * std::sqrt(2.) * sqr(r) * std::sqrt(n_I / (M_PI * chi_L_1p)) * sqr(1. + z) / std::sqrt(1. - z) / std::pow(denominator(z), 4);
    }
    static double series(std::vector<double> const& a, double z) {
        auto sum = 0.;
        for (auto n = 0u; n < a.size(); ++n) sum += a[n] * std::pow(z, n);
        return sum;
    }
    double g(double q2) const {
        auto zq = z(q2);
        return series(a_g, zq) / blaschke(zq, poles_1m) / phi_g(zq);
    }
    double f(double q2) const {
        auto zq = z(q2);
        return series(a_f, zq) / blaschke(zq, poles_1p) / phi_f(zq);
    }
    double F_1(double q2) const {
        auto zq = z(q2);
        return series(a_F1, zq) / blaschke(zq, poles_1p) / phi_F1(zq);
    }
    double F_2(double q2) const {
        auto zq = z(q2);
        return series(a_F2, zq) / blaschke(zq, poles_0m) / phi_F2(zq);
    }
    std::array<double, 4> operator()(double q2) const {
        auto a1 = f(q2) / (m_B + m_V);
        // F_1 = ((m_B + m_V) (m_B^2 - m_V^2 - q^2) A_1 - lambda A_2 / (m_B + m_V)) / (2 m_V)
        auto lambda = (sqr(m_B + m_V) - q2) * (sqr(m_B - m_V) - q2);
        auto a2 = (m_B + m_V) * ((m_B + m_V) * (sqr(m_B) - sqr(m_V) - q2) * a1 - 2. * m_V * F_1(q2)) / lambda;
        return {g(q2) * (m_B + m_V) / 2., F_2(q2) / 2., a1, a2};
    }
};

BGL const& B_D_star() {
    static BGL const bgl;
    return bgl;
}

// --- B_s -> D_s* and B_c -> J/psi: HPQCD z expansions of V, A_0, A_1, A_2,
//     F = 1 / P(q^2) sum_n a_n z(q^2, t_0)^n with P = prod_poles z(q^2, M_pole^2), t_+ = (m_B + m_D*)^2,
//     t_0 = (m_from - m_to)^2 and the B_c poles of 0-, 1-, 1+ states below the B D* threshold.
struct HPQCD {
    double t_plus, t_0;
    std::array<double, 4> V, A0, A1, A2;
    double F(std::array<double, 4> const& a, std::vector<double> const& poles, double q2) const {
        auto P = 1.;
        for (auto mass : poles) P *= z_map(q2, t_plus, sqr(mass));
        auto zq = z_map(q2, t_plus, t_0);
        auto sum = 0.;
        for (auto n = 0u; n < a.size(); ++n) sum += a[n] * std::pow(zq, n);
        return sum / P;
    }
    std::array<double, 4> operator()(double q2) const {
        static std::vector<double> const poles_0m{6.275, 6.872, 7.25}, poles_1m{6.335, 6.926, 7.02, 7.28}, poles_1p{6.745, 6.75, 7.15, 7.15};
        return {F(V, poles_1m, q2), F(A0, poles_0m, q2), F(A1, poles_1p, q2), F(A2, poles_1p, q2)};
    }
};

// B_s -> D_s*: HPQCD 21B (arXiv:2105.11433), eqs. (21)-(24) and tables XI-XIII. HPQCD 23 (arXiv:2304.03137) updates
// it, but publishes no coefficients for B_s -> D_s*; FLAG 2024 quotes only its F(1), which agrees.
HPQCD const& B_s_D_s_star() {
    static HPQCD const hpqcd{sqr(5.27964 + 2.010), sqr(5.3669 - 2.112),
        {0.100, -0.18, -0.006, 0.0}, {0.1047, -0.43, -0.10, -0.03}, {0.0552, -0.010, -0.03, 0.06}, {0.059, -0.11, -0.25, -0.05}};
    return hpqcd;
}

// B_c -> J/psi: HPQCD 20 (arXiv:2007.06957), eqs. (25)-(27), tables XIV-XVI, t_0 = t_-. HPQCD 25 (arXiv:2503.15090)
// supersedes it with a dispersive parametrisation whose coefficients are only in machine-readable supplementary files;
// its total rate agrees with HPQCD 20 within uncertainties.
HPQCD const& B_c_J_psi() {
    static HPQCD const hpqcd{sqr(5.27964 + 2.010), sqr(6.2749 - 3.0969),
        {0.1057, -0.746, 0.10, 0.006}, {0.1006, -0.731, 0.30, -0.02}, {0.0553, -0.266, 0.31, 0.11}, {0.0511, -0.22, -0.36, -0.05}};
    return hpqcd;
}

}

std::array<double, 2> pseudoscalar_form_factors(int id_from, int id_to, double q2, double m_from, double m_to) {
    id_from = std::abs(id_from);
    id_to = std::abs(id_to);
    auto bcl = [q2](BCL const& fit) -> std::array<double, 2> {
        return {fit.f_plus(q2), fit.f_0(q2)};
    };
    if (is_K(id_from) && is_pi(id_to)) return K_pi(id_from, q2);
    if (is_D(id_from) && is_K(id_to)) return bcl(D_K());
    if (is_D(id_from) && is_pi(id_to)) return bcl(D_pi());
    if (is_D_s(id_from) && is_eta(id_to)) return D_s_eta(q2);
    if (is_B(id_from) && is_pi(id_to)) return bcl(B_pi());
    if (is_B_s(id_from) && is_K(id_to)) return bcl(B_s_K());
    if (is_B(id_from) && is_D(id_to)) return bcl(B_D());
    if (is_B_s(id_from) && is_D_s(id_to)) return bcl(B_s_D_s());
    if (is_B_c(id_from) && id_to == 441) return B_c_eta_c(q2, m_from, m_to);
    print("no form factors for", id_from, "->", id_to);
    return {0., 0.};
}

namespace {

std::array<double, 4> fit_vector_form_factors(int id_from, int id_to, double q2) {
    if (is_D(id_from) && is_K_star(id_to)) return D_K_star(q2);
    if (is_B(id_from) && is_D_star(id_to)) return B_D_star()(q2);
    if (is_B(id_from) && is_rho(id_to)) return B_rho()(q2);
    if (is_B_s(id_from) && is_D_s_star(id_to)) return B_s_D_s_star()(q2);
    if (is_B_s(id_from) && is_K_star(id_to)) return B_s_K_star()(q2);
    if (is_B_c(id_from) && id_to == 443) return B_c_J_psi()(q2);
    print("no form factors for", id_from, "->", id_to);
    return {0., 0., 0., 0.};
}

}

double A_0_shift(int id_from, int id_to, double m_from, double m_to) {
    auto ff = fit_vector_form_factors(std::abs(id_from), std::abs(id_to), 0.);
    auto A_3 = ((m_from + m_to) * ff[2] - (m_from - m_to) * ff[3]) / (2. * m_to);
    return A_3 - ff[1];
}

std::array<double, 4> vector_form_factors(int id_from, int id_to, double q2, double m_from, double m_to) {
    auto ff = fit_vector_form_factors(std::abs(id_from), std::abs(id_to), q2);
    ff[1] += A_0_shift(id_from, id_to, m_from, m_to);
    return ff;
}

}
