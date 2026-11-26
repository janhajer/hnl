#pragma once

#include "Pythia8/Pythia.h"

namespace hnl {

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
}

}
