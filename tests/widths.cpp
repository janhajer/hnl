// Widths of MesonResonance and BaryonResonance for an HNL of negligible mass with |U|^2 = 1, which are the SM widths with a
// massless neutrino, against measured branching ratios (PDG 2025) and published lattice QCD rates.

#include <memory>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "BaryonDecays.hh"
#include "ResonanceWidths.hh"
#include "pythia.hh"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using hnl::sqr;

namespace {

double const hbar = 6.582119569e-25; // GeV s

struct Result {
    double width;    // GeV
    double V2;       // |V_CKM|^2 of Pythia, as set by hnl::set_physics_inputs
    double tau0;     // c tau of the parent [mm] before the HNL channels, as set by hnl::set_physics_inputs
};

// Gamma(parent -> [daughter] lepton N); daughter 0: two-body decay
Result hnl_width(int parent, int daughter, int lepton, int up, int down) {
    Pythia8::Pythia pythia(hnl::xml_dir(), false);
    hnl::set_pythia_branching_ratios(pythia);
    for (auto hnl : hnl::heavy_neutral_leptons()) hnl::set_pythia_stable(pythia, hnl, 1e-4);
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
        return {channel.bRatio() * particle.mWidth(), pythia.coupSM.V2CKMgen(up, down), tau0};
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
