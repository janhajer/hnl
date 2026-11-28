#pragma once

#include <algorithm>
#include <utility>
#include <vector>

#include "Pythia8/Pythia.h"

namespace hnl {

namespace detail {

inline bool has_products(Pythia8::DecayChannel& channel, std::vector<int> products) {
    if (channel.multiplicity() != static_cast<int>(products.size())) return false;
    std::vector<int> present;
    for (int i = 0; i < channel.multiplicity(); ++i) present.push_back(channel.product(i));
    std::sort(present.begin(), present.end());
    std::sort(products.begin(), products.end());
    return present == products;
}

inline bool is_semitauonic(Pythia8::DecayChannel& channel) { // tau+ nu_tau X
    bool tau = false, neutrino = false;
    for (int i = 0; i < channel.multiplicity(); ++i) {
        tau = tau || channel.product(i) == -15;
        neutrino = neutrino || channel.product(i) == 16;
    }
    return tau && neutrino;
}

// Sets the listed channels of a particle to measured branching ratios, scales its other tau+ nu_tau X channels so that
// all of them add up to the measured inclusive one (if given), and the remaining channels so that everything adds up
// to one. Running it again changes nothing.
inline void set_tau_branching_ratios(Pythia8::ParticleData& particle_data, int id,
                                     std::vector<std::pair<std::vector<int>, double>> const& measured, double inclusive = 0.) {
    auto entry = particle_data.particleDataEntryPtr(id);
    enum class Kind {listed, tau, other};
    std::vector<Kind> kinds;
    double listed = 0., tau = 0., other = 0.;
    for (int i = 0; i < entry->sizeChannels(); ++i) {
        auto& channel = entry->channel(i);
        auto found = std::find_if(measured.begin(), measured.end(), [&](auto const& m) { return has_products(channel, m.first); });
        if (found != measured.end()) channel.bRatio(found->second);
        kinds.push_back(found != measured.end() ? Kind::listed : is_semitauonic(channel) ? Kind::tau : Kind::other);
        (kinds.back() == Kind::listed ? listed : kinds.back() == Kind::tau ? tau : other) += channel.bRatio();
    }
    double tau_scale = inclusive > 0. && tau > 0. ? std::max(0., inclusive - listed) / tau : 1.;
    double other_scale = (1. - listed - tau * tau_scale) / other;
    for (int i = 0; i < entry->sizeChannels(); ++i) {
        if (kinds[i] == Kind::tau) entry->channel(i).rescaleBR(tau_scale);
        if (kinds[i] == Kind::other) entry->channel(i).rescaleBR(other_scale);
    }
}

}

// Inputs that replace Pythia 8.318 defaults. The resonances call this in their constructors, before Pythia::init() and
// before they read the SM lifetimes; it may run several times.
inline void set_physics_inputs(Pythia8::Pythia& pythia) {
    // |V_CKM| of the PDG 2025 global fit (CKM review, eq. (12.27)). Pythia's defaults are an older fit with larger
    // V_ub (0.00396) and V_cb (0.04221). Pythia limits V_cb to [0.0418, 0.0426], so the exclusive value cannot be used.
    pythia.readString("StandardModel:Vud = 0.97435");
    pythia.readString("StandardModel:Vus = 0.22501");
    pythia.readString("StandardModel:Vub = 0.003732");
    pythia.readString("StandardModel:Vcd = 0.22487");
    pythia.readString("StandardModel:Vcs = 0.97349");
    pythia.readString("StandardModel:Vcb = 0.04183");
    pythia.readString("StandardModel:Vtd = 0.00858");
    pythia.readString("StandardModel:Vts = 0.04111");
    pythia.readString("StandardModel:Vtb = 0.999118");
    // c tau [mm] of PDG 2025 for the parents of HNLs: the HNL branching ratios are Gamma(X -> N) times these lifetimes.
    // Pythia differs by -9.7 % (B_c), -3.4 % (B_s), +0.9 % (B0) and +0.7 % (D+); the others agree to 0.3 %.
    auto& particle_data = pythia.particleData;
    particle_data.tau0(211, 7804.5);  // pi+, 2.6033e-8 s
    particle_data.tau0(321, 3711.4);  // K+, 1.2380e-8 s
    particle_data.tau0(130, 15337.);  // K_L, 5.116e-8 s
    particle_data.tau0(310, 26.844);  // K_S, 0.8954e-10 s
    particle_data.tau0(411, 0.30969); // D+, 1033 fs
    particle_data.tau0(421, 0.12300); // D0, 410.3 fs
    particle_data.tau0(431, 0.15026); // D_s+, 501.2 fs
    particle_data.tau0(511, 0.45479); // B0, 1.517 ps
    particle_data.tau0(521, 0.49106); // B+, 1.638 ps
    particle_data.tau0(531, 0.45449); // B_s, 1.516 ps
    particle_data.tau0(541, 0.15289); // B_c+, 0.510 ps
    particle_data.tau0(5122, 0.44010); // Lambda_b, 1.468 ps
    // Branching ratios into tau of PDG 2025, for HNLs from tau decays. Pythia has D_s -> tau nu 6.40 % and B+, B0 -> tau nu X
    // 3.0 %, of which 0.66 % are the unmeasured D** tau nu. D_s -> tau nu 5.39(9) %, D+ -> tau nu 1.20(27)e-3, B+ -> tau nu
    // 1.09(24)e-4; B -> D tau nu 0.83(8) % and D* tau nu 1.42(7) % (B+/B0 admixture) and b-hadron -> tau nu X 2.41(23) %
    // (b admixture, mean life 1.5673 ps), shared between B+ and B0 by their lifetimes; the other tau nu X channels are
    // scaled to the inclusive rate (by about 0.25). B_s, B_c and Lambda_b keep Pythia's values (no measurements).
    detail::set_tau_branching_ratios(particle_data, 431, {{{-15, 16}, 0.0539}});
    detail::set_tau_branching_ratios(particle_data, 411, {{{-15, 16}, 0.00120}});
    detail::set_tau_branching_ratios(particle_data, 521, {{{16, -15, -421}, 0.00862}, {{16, -15, -423}, 0.01474}, {{16, -15}, 1.09e-4}}, 0.02519);
    detail::set_tau_branching_ratios(particle_data, 511, {{{16, -15, -411}, 0.00798}, {{16, -15, -413}, 0.01366}}, 0.02333);
}

}
