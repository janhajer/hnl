#pragma once

#include <array>
#include <functional>

#include "Pythia8/ResonanceWidths.h"

namespace Pythia8{
    class Pythia;
}

namespace hnl
{

// Lambda_b -> Lambda_c l N with the lattice QCD form factors of Detmold, Lehner, Meinel, arXiv:1503.01421 (nominal fit,
// helicity basis, physical limit). The spin-summed |M|^2 comes from Dirac traces of the hadronic and leptonic currents,
// so it holds for any masses of l and N.

// {f_+, f_0, f_perp, g_+, g_0, g_perp} at q2 [GeV^2]
std::array<double, 6> lambda_b_form_factors(double q2);

// sum over spins of |M|^2 / (G_F^2 |V_cb|^2 |U|^2), averaged over the Lambda_b spin, for
// Lambda_b(p) -> Lambda_c(p_c) l-(k_l) N(k_N), M = G_F / sqrt(2) V_cb U <Lambda_c|J^mu|Lambda_b> [u_l gamma_mu (1 - gamma_5) v_N]
double lambda_b_matrix_element(Pythia8::Vec4 const& p, Pythia8::Vec4 const& p_c, Pythia8::Vec4 const& k_l, Pythia8::Vec4 const& k_N);

// Gamma(Lambda_b -> Lambda_c l N) / (G_F^2 |V_cb|^2 |U|^2) [GeV^5]
double lambda_b_width(double m_b, double m_c, double m_l, double m_N);

// Lambda_b with its HNL channels Lambda_b -> Lambda_c+ l- N; width SM plus HNL, as MesonResonance
struct BaryonResonance : public Pythia8::ResonanceWidths {
    BaryonResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int id_from = 5122);
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
