// HNL production in tau decays (TauResonance): the three-body widths against the integrals of 1805.08567 and the massless
// limit, the two-body widths of a nearly massless HNL against the measured BR(tau -> nu pi) and BR(tau -> nu K) (PDG 2025).

#include <algorithm>
#include <cmath>
#include <functional>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "Pythia8/Pythia.h"

#include "PhysicsInputs.hh"
#include "TauResonance.hh"
#include "pythia.hh"

using Catch::Matchers::WithinRel;

namespace {

double const m_tau = 1.77686, m_mu = 0.1056584;

double kallen(double a, double b, double c) {
    return a * a + b * b + c * c - 2. * a * b - 2. * a * c - 2. * b * c;
}

double integrate(std::function<double(double)> const& f, double low, double high) {
    return boost::math::quadrature::gauss_kronrod<double, 61>::integrate(f, low, high, 15, 1e-12);
}

// Gamma(tau -> N l nubar_l) / (G_F^2 |U_tau|^2) of 1805.08567
double paper_N_l_nu(double m_N, double m_l) {
    double y_N = m_N / m_tau, y_l = m_l / m_tau;
    auto f = [&](double x) {
        return std::pow(x - y_l * y_l, 2) * std::sqrt(std::max(0., kallen(1., x, y_N * y_N))) / std::pow(x, 3)
            * ((x + 2. * y_l * y_l) * std::pow(1. - y_N * y_N, 2) + x * (x - y_l * y_l) * (1. + y_N * y_N - y_l * y_l) - x * std::pow(y_l, 4) - 2. * std::pow(x, 3));
    };
    return std::pow(m_tau, 5) / (96. * std::pow(M_PI, 3)) * integrate(f, y_l * y_l, std::pow(1. - y_N, 2));
}

// Gamma(tau -> nu_tau l N) / (G_F^2 |U_l|^2) of 1805.08567
double paper_nu_l_N(double m_N, double m_l) {
    double y_N = m_N / m_tau, y_l = m_l / m_tau;
    auto f = [&](double x) {
        return std::pow(1. - x, 2) * std::sqrt(std::max(0., kallen(x, y_N * y_N, y_l * y_l))) / std::pow(x, 3)
            * (2. * std::pow(x, 3) + x - x * (1. - x) * (1. - y_N * y_N - y_l * y_l) - (2. + x) * std::pow(y_N * y_N - y_l * y_l, 2));
    };
    return std::pow(m_tau, 5) / (96. * std::pow(M_PI, 3)) * integrate(f, std::pow(y_l + y_N, 2), 1.);
}

}

TEST_CASE("tau -> l nu nu with massless products gives m_tau^5 / (192 pi^3)", "[tau]") {
    CHECK_THAT(hnl::tau_leptonic_width(m_tau, 0., 0., 0.), WithinRel(std::pow(m_tau, 5) / (192. * std::pow(M_PI, 3)), 1e-8));
    CHECK(hnl::tau_leptonic_width(m_tau, 0., m_mu, 500.) == 0.); // closed channel: an HNL heavier than the tau
    CHECK(hnl::tau_leptonic_width(m_tau, 1.7, m_mu, 0.) == 0.);
}

TEST_CASE("tau -> N l nubar and tau -> nu_tau l N reproduce the integrals of 1805.08567", "[tau]") {
    for (double m_N : {0.1, 0.5, 1.0, 1.5}) {
        CHECK_THAT(hnl::tau_leptonic_width(m_tau, m_N, m_mu, 0.), WithinRel(paper_N_l_nu(m_N, m_mu), 1e-6));
        CHECK_THAT(hnl::tau_leptonic_width(m_tau, 0., m_mu, m_N), WithinRel(paper_nu_l_N(m_N, m_mu), 1e-6));
    }
}

TEST_CASE("tau -> N pi and N K of a light HNL reproduce the measured tau -> nu pi and nu K", "[tau]") {
    Pythia8::Pythia pythia(hnl::xml_dir(), false);
    hnl::set_physics_inputs(pythia);
    pythia.readString("Print:quiet = on");
    pythia.init();
    double width = hnl::tau_to_Gamma(pythia.particleData.tau0(15));
    double G_F2 = std::pow(pythia.coupSM.GF(), 2);
    double pi = G_F2 * pythia.coupSM.V2CKMgen(1, 1) * hnl::tau_meson_width(m_tau, 1e-4, 211, pythia.particleData.m0(211));
    double K = G_F2 * pythia.coupSM.V2CKMgen(1, 2) * hnl::tau_meson_width(m_tau, 1e-4, 321, pythia.particleData.m0(321));
    CHECK_THAT(pi / width, WithinRel(0.1082, 0.02)); // PDG 2025: 10.82(5)%
    CHECK_THAT(K / width, WithinRel(0.00696, 0.03)); // PDG 2025: 0.696(10)%
}

TEST_CASE("tau branching ratios of D_s, D+, B+ and B0 are those of PDG 2025", "[tau]") {
    Pythia8::Pythia pythia(hnl::xml_dir(), false);
    hnl::set_physics_inputs(pythia);
    hnl::set_physics_inputs(pythia); // the resonance constructors call it several times
    auto sum = [&](int id, bool semitauonic_only) {
        double result = 0.;
        auto entry = pythia.particleData.particleDataEntryPtr(id);
        for (int i = 0; i < entry->sizeChannels(); ++i) {
            auto const& channel = entry->channel(i);
            bool tau = false, neutrino = false;
            for (int j = 0; j < channel.multiplicity(); ++j) {
                tau = tau || channel.product(j) == -15;
                neutrino = neutrino || channel.product(j) == 16;
            }
            if (!semitauonic_only || (tau && neutrino)) result += channel.bRatio();
        }
        return result;
    };
    CHECK_THAT(sum(431, true), WithinRel(0.0539, 1e-9));
    CHECK_THAT(sum(411, true), WithinRel(0.00120, 1e-9));
    CHECK_THAT(sum(521, true), WithinRel(0.02519, 1e-9));
    CHECK_THAT(sum(511, true), WithinRel(0.02333, 1e-9));
    for (int id : {431, 411, 521, 511}) CHECK_THAT(sum(id, false), WithinRel(1., 1e-9));
}
