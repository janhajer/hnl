#include "TauResonance.hh"

#include <algorithm>
#include <cmath>

#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "Pythia8/Pythia.h"

#include "PhysicsInputs.hh"
#include "ResonanceWidths.hh"
#include "id.hh"
#include "math.hh"
#include "pythia.hh"

namespace hnl {

namespace {

double kallen(double a, double b, double c) {
    return sqr(a) + sqr(b) + sqr(c) - 2. * a * b - 2. * a * c - 2. * b * c;
}

bool is_vector_meson(int meson) {
    return meson == 213 || meson == 323;
}

}

// Gamma(tau -> N h_P) and Gamma(tau -> N h_V) of 1805.08567 without G_F^2 |V_UD|^2 |U_tau|^2
double tau_meson_width(double m_tau, double m_N, int meson, double m_meson) {
    if (m_tau <= m_N + m_meson) return 0.;
    double y_N2 = sqr(m_N / m_tau), y_h2 = sqr(m_meson / m_tau);
    double root = std::sqrt(kallen(1., y_N2, y_h2));
    double constant = decay_constant(meson);
    if (is_vector_meson(meson))
        return sqr(constant) * Pythia8::pow3(m_tau) / (16. * M_PI * sqr(m_meson)) * (sqr(1. - y_N2) + y_h2 * (1. + y_N2 - 2. * y_h2)) * root;
    return sqr(constant) * Pythia8::pow3(m_tau) / (16. * M_PI) * (sqr(1. - y_N2) - y_h2 * (1. + y_N2)) * root;
}

double tau_leptonic_matrix_element(Pythia8::Vec4 const& p, Pythia8::Vec4 const& k_a, Pythia8::Vec4 const& k_b, Pythia8::Vec4 const& k_c) {
    return 64. * (p * k_c) * (k_a * k_b);
}

// dGamma / dm_ab^2 = |M|^2(m_ab^2) (m_bc^2 max - m_bc^2 min) / ((2 pi)^3 32 m_tau^3); |M|^2 depends on m_ab^2 only
double tau_leptonic_width(double m_tau, double m_a, double m_b, double m_c) {
    if (m_tau <= m_a + m_b + m_c) return 0.;
    double low = sqr(m_a + m_b), high = sqr(m_tau - m_c);
    auto integrand = [&](double m_ab2) {
        double m_ab = std::sqrt(m_ab2);
        double e_b = (m_ab2 - sqr(m_a) + sqr(m_b)) / (2. * m_ab); // in the a b rest frame
        double e_c = (sqr(m_tau) - m_ab2 - sqr(m_c)) / (2. * m_ab);
        double range = 4. * std::sqrt(std::max(0., sqr(e_b) - sqr(m_b))) * std::sqrt(std::max(0., sqr(e_c) - sqr(m_c)));
        double p_k_c = (sqr(m_tau) + sqr(m_c) - m_ab2) / 2.;
        double k_a_k_b = (m_ab2 - sqr(m_a) - sqr(m_b)) / 2.;
        return 64. * p_k_c * k_a_k_b * range;
    };
    double integral = boost::math::quadrature::gauss_kronrod<double, 61>::integrate(integrand, low, high, 15, 1e-10);
    return integral / (Pythia8::pow3(2. * M_PI) * 32. * Pythia8::pow3(m_tau));
}

TauResonance::TauResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_) :
    neutrino_coupling(neutrino_coupling_) {
    set_physics_inputs(pythia);
    initBasic(15);
    particlePtr = pythia.particleData.particleDataEntryPtr(idRes);
    tau0_sm = particlePtr->tau0();
}

bool TauResonance::init(Pythia8::Info* infoPtrIn) {
    sum = 0.;
    if (!ResonanceWidths::init(infoPtrIn)) return false;
    // as MesonResonance: width SM plus HNL channels, the channels open for the particle add up to one
    double width = tau_to_Gamma(tau0_sm) + sum;
    particlePtr->rescaleBR(particlePtr->mWidth() / width);
    particlePtr->setMWidth(width, false);
    particlePtr->setTau0(Gamma_to_tau(width));
    return true;
}

bool TauResonance::initBSM() {
    particlePtr->setTau0(tau0_sm); // the previous init() has changed it
    particlePtr->setMWidth(tau_to_Gamma(tau0_sm));
    particlePtr->setIsResonance(false);
    particlePtr->setMayDecay(true);
    return true;
}

// tau- -> N pi-, N K-, N rho-, N K*-, N l- nubar_l for the HNL mixing with nu_tau, tau- -> nu_tau l- N for those mixing
// with nu_e or nu_mu; channels of HNLs without that mixing get zero width
bool TauResonance::allowCalc() {
    if (channels_added) return true; // init() runs in the constructor of Pythia's resonance list and in pythia.init()
    for (auto meson : {211, 321, 213, 323}) particlePtr->addChannel(1, 0., 0, heavy_neutral_tau, -meson);
    for (auto lepton : {11, 13}) particlePtr->addChannel(1, 0., 0, heavy_neutral_tau, lepton, -(lepton + 1));
    particlePtr->addChannel(1, 0., 0, 16, 11, heavy_neutral_electron);
    particlePtr->addChannel(1, 0., 0, 16, 13, heavy_neutral_muon);
    channels_added = true;
    return true;
}

void TauResonance::calcWidth(bool) {
    auto const& channel = particlePtr->channel(iChannel);
    if (!has_neutrino(channel)) { // SM channel, keeps its width
        widNow = GammaRes * channel.bRatio();
        return;
    }
    widNow = 0.;
    double m_tau = particleDataPtr->m0(idRes);
    double G_F2 = sqr(coupSMPtr->GF());
    if (channel.multiplicity() == 2) { // N h
        int neutrino = channel.product(0), meson = std::abs(channel.product(1));
        double coupling = neutrino_coupling(neutrino, 16);
        if (coupling <= 0.) return;
        double V2 = meson == 211 || meson == 213 ? coupSMPtr->V2CKMgen(1, 1) : coupSMPtr->V2CKMgen(1, 2);
        widNow = G_F2 * V2 * coupling * tau_meson_width(m_tau, particleDataPtr->m0(neutrino), meson, particleDataPtr->m0(meson));
    } else if (is_heavy_neutral_lepton(channel.product(0))) { // N l nubar_l
        int neutrino = channel.product(0), lepton = channel.product(1);
        double coupling = neutrino_coupling(neutrino, 16);
        if (coupling <= 0.) return;
        widNow = G_F2 * coupling * tau_leptonic_width(m_tau, particleDataPtr->m0(neutrino), particleDataPtr->m0(lepton), 0.);
    } else { // nu_tau l N
        int lepton = channel.product(1), neutrino = channel.product(2);
        double coupling = neutrino_coupling(neutrino, lepton + 1);
        if (coupling <= 0.) return;
        widNow = G_F2 * coupling * tau_leptonic_width(m_tau, 0., particleDataPtr->m0(lepton), particleDataPtr->m0(neutrino));
    }
    sum += widNow;
}

}
