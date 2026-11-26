# Heavy Neutral Leptons in Pythia8

## Dependencies

* Pythia 8.3
* Boost
    * Filesystem
    * Range
    * Optional
    * Algorithm
    * Io

## HNLs from B mesons

`b-mesons <m_N> <|U_mu|^2> <events> <mode>` (`main-b-mesons.cpp`; mode 0: SM branching ratios, 1: all B into HNLs, 2: one B per event into HNLs) shows the setup:

* `NeutrinoResonance`: HNL mass, widths and decays (Majorana by default),
* `MesonResonance` for each meson: HNL channels added to its decay table, widths from arXiv:1805.08567,
* `MesonDecays`: decay handler that generates the HNL channels with the matrix element instead of phase space.
* `BaryonResonance`: Lambda_b -> Lambda_c l N with the lattice QCD form factors of arXiv:1503.01421; |M|^2 from Dirac traces, for any lepton and HNL masses. In the SM limit it gives their Gamma / |V_cb|^2 for e, mu and tau (21.5, 21.5, 7.15 ps^-1).

## HNL decays

Widths of arXiv:1805.08567. Hadronic decays:

* one meson at all masses: `l P` (pi, K, D, D_s, B, B_c), `l V` (rho, K\*, D\*, D_s\*, a_1), `nu P` (pi0, eta, eta', eta_c), `nu V` (rho0, omega, phi, J/psi, a_1). The K\* and a_1 couplings come from tau decays, as g_rho in the paper. The neutral-current factors of omega, phi and J/psi follow from the quark couplings with the paper's decay constants and are smaller than in its table (see `correction_factor`);
* for each quark current, the remainder of the inclusive width (quarks with the QCD correction, alpha_s at three loops) over its one-meson channels, if positive and above the lightest states those channels do not cover. The quark pair fragments in Pythia's hadron level (meMode 94 with the V-A matrix element for the charged current, 93 for the neutral one).

At m_N = m_tau the charged-current channels reproduce tau decays (pi : rho : a_1 : rest).

Decay HNLs with `next()` or `forceHadronLevel()`. `moreDecays()` with its default `allowPartons = false` replaces the quark pair by a random-flavour pair: the charge is not conserved, and below 1.45 GeV `ParticleDecays::pickHadrons` loops forever.

`patches/pythia8318-ParticleDecays.patch` fixes this in Pythia 8.318 (`patch -p1 < .../pythia8318-ParticleDecays.patch` in the Pythia directory, then `make`): `moreDecays(false)` keeps the quark flavours and the flavour loop has a try limit; `next()` is not affected. HNL decays at 0.5-4.5 GeV then conserve the charge and never hang; so do SM decays, where unpatched `moreDecays()` gives the wrong charge in 43% of B+ and 23% of Lambda_b decays. Up to m_N = 2 GeV the multiplicity agrees with `forceHadronLevel()` to 2%; above, `moreDecays()` picks the hadrons by phase space (charged multiplicity +22% at 3 GeV), so `forceHadronLevel()` stays the better choice there.
