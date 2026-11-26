// Widths of MesonResonance and BaryonResonance for an HNL of negligible mass with |U|^2 = 1, which are the SM widths with a
// massless neutrino, against measured branching ratios (PDG 2025) and published lattice QCD rates; and the B -> D(*) l N
// widths of a heavy HNL against their zero-recoil limit.

#include <cmath>
#include <memory>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "BaryonDecays.hh"
#include "FormFactors.hh"
#include "ResonanceWidths.hh"
#include "pythia.hh"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using hnl::cube;
using hnl::sqr;

namespace {

double const hbar = 6.582119569e-25; // GeV s

struct Result {
    double width;    // GeV
    double V2;       // |V_CKM|^2 of Pythia, as set by hnl::set_physics_inputs
    double tau0;     // c tau of the parent [mm] before the HNL channels, as set by hnl::set_physics_inputs
    double GF;       // GeV^-2
    double m_parent; // masses [GeV] of Pythia's particle data, used by the code
    double m_daughter;
    double m_lepton;
};

// Gamma(parent -> [daughter] lepton N) for an HNL of mass m_N [GeV]; daughter 0: two-body decay
Result hnl_width(int parent, int daughter, int lepton, int up, int down, double m_N = 1e-4) {
    Pythia8::Pythia pythia(hnl::xml_dir(), false);
    hnl::set_pythia_branching_ratios(pythia);
    for (auto hnl : hnl::heavy_neutral_leptons()) hnl::set_pythia_stable(pythia, hnl, m_N);
    auto resonance = std::make_shared<hnl::MesonResonance>(pythia, hnl::neutrino_coupling(1.), parent);
    auto tau0 = pythia.particleData.tau0(parent);
    pythia.addResonancePtr(resonance);
    pythia.init();
    auto& particle = *pythia.particleData.particleDataEntryPtr(parent);
    auto hnl = lepton + 9900001;
    for (auto i = 0; i < particle.sizeChannels(); ++i) {
        auto const& channel = particle.channel(i);
        if (channel.multiplicity() != (daughter == 0 ? 2 : 3)) continue;
        if (std::abs(channel.product(0)) != hnl || std::abs(channel.product(1)) != lepton) continue;
        if (daughter != 0 && std::abs(channel.product(2)) != daughter) continue;
        return {channel.bRatio() * particle.mWidth(), pythia.coupSM.V2CKMgen(up, down), tau0, pythia.coupSM.GF(),
                particle.m0(), daughter == 0 ? 0. : pythia.particleData.m0(daughter), pythia.particleData.m0(lepton)};
    }
    FAIL("channel not found");
    return {};
}

// branching ratio with the lifetime of hnl::set_physics_inputs
double branching_ratio(Result const& result) {
    return result.width / hnl::tau_to_Gamma(result.tau0);
}

// Gamma / |V|^2 [1/s]
double reduced_width(Result const& result) {
    return result.width / result.V2 / hbar;
}

double kallen(double a, double b, double c) {
    return sqr(a) + sqr(b) + sqr(c) - 2 * (a * b + a * c + b * c);
}

// Zero-recoil limit of Gamma(B -> H l N) for the SM current with |U|^2 = 1 (unpublished calculation by S. Biondini and
// J. Hajer): with w = v.v', for a single HNL mass eigenstate (Dirac or Majorana)
//   dGamma / dw = m_H^3 / (2 pi^2) sqrt(w^2 - 1) rho(q^2),
//   rho_F  = lambda^1/2 / (4 pi q^4) [(q^2 - m_+^2) |m_- a_V|^2 + (q^2 - m_-^2) |m_+ b_V|^2],
//   rho_GT = lambda^1/2 / (4 pi q^4) [(q^2 - m_-^2) (2 q^2 + m_+^2) |a_A|^2 + (q^2 - m_+^2) (2 q^2 + m_-^2) |b_A|^2],
// m_+- = m_l +- m_N, lambda = lambda(q^2, m_l^2, m_N^2). The Fermi (vector) source produces only H = D, the Gamow-Teller
// (axial-vector) source only H = D*. For the SM current J_V = -J_A = G_F V U / sqrt2 l gamma (1 - gamma_5) N, so
// a_V = -b_V = -a_A = b_A = eta G_F V U / sqrt2, with the zero-recoil normalisations taken from the form factors of the
// code at q^2_max = (m_B - m_H)^2: eta_V = (m_B + m_D) f_0 / (2 sqrt(m_B m_D)) and eta_A = (m_B + m_D*) A_1 / (2 sqrt(m_B
// m_D*)). The limit holds up to relative corrections of order w_max - 1.
double zero_recoil_width(Result const& r, int parent, int daughter, double m_N) {
    auto m_B = r.m_parent;
    auto m_H = r.m_daughter;
    auto m_l = r.m_lepton;
    auto q2_max = sqr(m_B - m_H);
    auto vector = hnl::is_vector(daughter);
    auto eta = (m_B + m_H) / 2 / std::sqrt(m_B * m_H) * (vector
        ? hnl::vector_form_factors(parent, daughter, q2_max, m_B, m_H)[2]
        : hnl::pseudoscalar_form_factors(parent, daughter, q2_max, m_B, m_H)[1]);
    auto coupling2 = sqr(eta * r.GF) * r.V2 / 2;
    auto m_p2 = sqr(m_l + m_N);
    auto m_m2 = sqr(m_l - m_N);
    auto rho = [&](double q2) {
        auto bracket = vector
            ? (q2 - m_m2) * (2 * q2 + m_p2) + (q2 - m_p2) * (2 * q2 + m_m2)
            : (q2 - m_p2) * m_m2 + (q2 - m_m2) * m_p2;
        return std::sqrt(std::max(kallen(q2, sqr(m_l), sqr(m_N)), 0.)) / (4 * M_PI * sqr(q2)) * coupling2 * bracket;
    };
    auto w_max = (sqr(m_B) + sqr(m_H) - m_p2) / (2 * m_B * m_H);
    // w = 1 + (w_max - 1) (1 - cos theta) / 2 removes the square-root end points
    auto n = 2000;
    auto sum = 0.;
    for (auto i = 0; i < n; ++i) {
        auto theta = M_PI * (i + .5) / n;
        auto w = 1 + (w_max - 1) * (1 - std::cos(theta)) / 2;
        auto q2 = sqr(m_B) + sqr(m_H) - 2 * m_B * m_H * w;
        sum += (w_max - 1) / 2 * std::sin(theta) * std::sqrt(sqr(w) - 1) * rho(q2);
    }
    return cube(m_H) / (2 * sqr(M_PI)) * sum * M_PI / n;
}

// HNL mass for which w_max = 1 + delta
double endpoint_mass(Result const& r, double delta) {
    auto q2_min = sqr(r.m_parent) + sqr(r.m_daughter) - 2 * r.m_parent * r.m_daughter * (1 + delta);
    return std::sqrt(q2_min) - r.m_lepton;
}

// Gamma / Gamma_zero-recoil - 1 at w_max = 1 + delta
double endpoint_deviation(int parent, int daughter, int lepton, double delta) {
    auto masses = hnl_width(parent, daughter, lepton, 2, 3);
    auto m_N = endpoint_mass(masses, delta);
    auto result = hnl_width(parent, daughter, lepton, 2, 3, m_N);
    return result.width / zero_recoil_width(result, parent, daughter, m_N) - 1;
}

}

TEST_CASE("Pythia inputs are replaced", "[inputs]") {
    auto result = hnl_width(541, 443, 13, 2, 3);
    CHECK_THAT(result.V2, WithinRel(sqr(0.04183), 1e-6));
    CHECK_THAT(result.tau0, WithinRel(0.15289, 1e-6));
}

// Measured branching ratios include radiative and short-distance electroweak corrections that the code does not have
// (about 1-3 %); the tolerances allow for these and for the uncertainties of the form factors.

TEST_CASE("two-body decays reproduce the measured branching ratios", "[widths]") {
    CHECK_THAT(branching_ratio(hnl_width(211, 0, 13, 1, 1)), WithinRel(0.99988, 3e-2));
    CHECK_THAT(branching_ratio(hnl_width(321, 0, 13, 1, 2)), WithinRel(0.6356, 3e-2));
    CHECK_THAT(branching_ratio(hnl_width(431, 0, 13, 2, 2)), WithinRel(5.35e-3, 5e-2));
}

TEST_CASE("K and D semileptonic decays reproduce the measured branching ratios", "[widths]") {
    CHECK_THAT(branching_ratio(hnl_width(130, 211, 11, 1, 2)), WithinRel(0.4055, 3e-2));
    CHECK_THAT(branching_ratio(hnl_width(321, 111, 11, 1, 2)), WithinRel(0.0507, 3e-2));
    CHECK_THAT(branching_ratio(hnl_width(421, 321, 11, 2, 2)), WithinRel(0.03538, 2e-2));
    CHECK_THAT(branching_ratio(hnl_width(421, 211, 11, 2, 1)), WithinRel(2.91e-3, 3e-2));
    // normalised to these two
    CHECK_THAT(branching_ratio(hnl_width(411, 313, 11, 2, 2)), WithinRel(0.0540, 1e-2));
    CHECK_THAT(branching_ratio(hnl_width(431, 221, 11, 2, 2)), WithinRel(0.0227, 1e-2));
}

TEST_CASE("B semileptonic decays reproduce the measured branching ratios with the |V| of their fits", "[widths]") {
    auto tau_B0 = 1.517e-12;
    // FLAG 2024, tables 50, 52 and 53: |V_ub| = 3.61e-3, |V_cb| = 40.0e-3 and 39.23e-3; PDG 2025 for l = e, mu
    CHECK_THAT(reduced_width(hnl_width(511, 211, 11, 1, 3)) * sqr(3.61e-3) * tau_B0, WithinRel(1.50e-4, 8e-2));
    CHECK_THAT(reduced_width(hnl_width(511, 411, 11, 2, 3)) * sqr(40.0e-3) * tau_B0, WithinRel(0.0210, 8e-2));
    CHECK_THAT(reduced_width(hnl_width(511, 413, 11, 2, 3)) * sqr(39.23e-3) * tau_B0, WithinRel(0.0487, 4e-2));
}

TEST_CASE("B_s -> D_s* and B_c -> J/psi reproduce HPQCD", "[widths]") {
    // Gamma / |eta_EW V_cb|^2: HPQCD 21B (arXiv:2105.11433) for l = e, HPQCD 20 (arXiv:2007.06957) for l = mu
    CHECK_THAT(reduced_width(hnl_width(531, 433, 11, 2, 3)), WithinRel(2.07e13, 3e-2));
    CHECK_THAT(reduced_width(hnl_width(541, 443, 13, 2, 3)), WithinRel(1.73e13, 3e-2));
}

TEST_CASE("Lambda_b -> Lambda_c l nu reproduces Detmold, Lehner, Meinel", "[widths]") {
    // Gamma / |V_cb|^2 = 21.5 (e, mu) and 7.15 (tau) ps^-1, arXiv:1503.01421
    auto width = [](double m_l) {
        return sqr(1.1663788e-5) * hnl::lambda_b_width(5.6196, 2.28646, m_l, 0.) / hbar * 1e-12;
    };
    CHECK_THAT(width(0.000511), WithinRel(21.5, 1e-2));
    CHECK_THAT(width(1.77686), WithinRel(7.15, 1e-2));
}

TEST_CASE("B -> D(*) l N approaches its zero-recoil limit at the kinematic end point", "[widths]") {
    // The deviation vanishes linearly in w_max - 1, with slopes of about 1.4 (D) and 0.6 (D*). Its Richardson extrapolation
    // 2 d(delta / 2) - d(delta) to w_max = 1, which removes the linear term, tests the normalisation.
    for (auto [daughter, lepton] : {std::pair{411, 13}, {411, 15}, {413, 13}, {413, 15}}) {
        INFO("B0 -> " << daughter << " " << lepton << " N");
        auto coarse = endpoint_deviation(511, daughter, lepton, 0.01);
        auto fine = endpoint_deviation(511, daughter, lepton, 0.005);
        CHECK(std::abs(coarse) < 2 * 0.01);
        CHECK(std::abs(fine) < 2 * 0.005);
        CHECK_THAT(2 * fine - coarse, WithinAbs(0., 1e-3));
    }
}
