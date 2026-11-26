#pragma once

#include <functional>

#include "Pythia8/Basics.h"
#include "Pythia8/ResonanceWidths.h"

namespace Pythia8 {
class Pythia;
}

namespace hnl {

// HNL production in tau decays, arXiv:1805.08567 (subsection 'HNL production from tau lepton'): tau -> N h (h = pi, K, rho, K*) and tau -> N l nubar_l
// for mixing with nu_tau, tau -> nu_tau l N for mixing with nu_e or nu_mu (l = e, mu).

// Gamma(tau -> N h) / (G_F^2 |V_UD|^2 |U_tau|^2) [GeV^3]; f_h for pseudoscalars, g_h = f_h m_h for vectors
double tau_meson_width(double m_tau, double m_N, int meson, double m_meson);

// spin-averaged |M|^2 / (G_F^2 |U|^2) = 64 (p.k_c)(k_a.k_b) of tau(p) -> a b c through a virtual W, with (a, b, c) =
// (N, l, nubar_l) for mixing with nu_tau and (nu_tau, l, N) for mixing with nu_l
double tau_leptonic_matrix_element(Pythia8::Vec4 const& p, Pythia8::Vec4 const& k_a, Pythia8::Vec4 const& k_b, Pythia8::Vec4 const& k_c);

// Gamma(tau -> a b c) / (G_F^2 |U|^2) [GeV^5] of that matrix element; m_tau^5 / (192 pi^3) for massless products
double tau_leptonic_width(double m_tau, double m_a, double m_b, double m_c);

// tau with its HNL channels; width SM plus HNL, as MesonResonance and BaryonResonance
struct TauResonance : public Pythia8::ResonanceWidths {
    TauResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_);
    bool init(Pythia8::Info* infoPtrIn) override;
protected:
    bool initBSM() override;
    bool allowCalc() override;
    void calcWidth(bool called_from_init = false) override;
private:
    std::function<double (int id_heavy, int id_light)> neutrino_coupling;
    Pythia8::ParticleDataEntryPtr particlePtr;
    double sum = 0.;
    bool channels_added = false;
    double tau0_sm;
};

}
