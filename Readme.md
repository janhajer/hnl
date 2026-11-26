# Heavy Neutral Leptons in Pythia8

## Dependencies

* Pythia 8.3 (tested with 8.318)
* Boost headers (Range, Optional, Algorithm), e.g. `libboost-dev`
* CMake >= 3.15
* HepMC3, optional, only for `write-hepmc`, e.g. `libhepmc3-dev`

## Build

```
mkdir build && cd build
cmake .. -DPYTHIA8_DIR=<pythia prefix>   # default: $PYTHIA8, else ~/Development/pythia
make -j
```

`write-hepmc` is built if CMake finds HepMC3; otherwise add `-DHepMC3_DIR=<prefix>/share/HepMC3/cmake`.

The executables read Pythia's data from `<pythia prefix>/share/Pythia8/xmldoc`, fixed at build time; the environment variable `PYTHIA8DATA` overrides it.

## HNLs in HepMC3

`write-hepmc <m_N> <|U_mu|^2> <flavour> [events] [seed]` (`write-hepmc.hh`) writes `<m_N>-<flavour>.hepmc` (HepMC3 ascii) with `events` (default 1000) events in which an HNL (`nu_Rmu`, mixing with nu_mu only, stable) from a parent of the given flavour has 1.4 < |eta| < 3.5:

| flavour | parents |
|---|---|
| `b` | B0, B+, Bs, Bc, Lambda_b, Upsilon |
| `c` | D0, D+, Ds, J/psi |
| `s` | K+, K_L, K_S |
| `pi` | pi+ |

Events are inclusive pp collisions (`SoftQCD:nonDiffractive`, as in `b-mesons`); as in mode 2 of `b-mesons`, one parent per event decays into an HNL, picked according to the branching ratios, and the event weight is the sum of these branching ratios. Only parents present after `next()` are found, so a parent must not come from the decay of another: the flavours are separate runs, and in the `c` run the b hadrons decay as in the SM and their D mesons are parents. The total is the sum over the flavours.

* Run info: the attributes `hnl_id`, `hnl_mass_GeV`, `hnl_U2` and `parent_flavour`.
* Events: the weight and the cross section [pb] of HNLs in the acceptance so far, sigma_gen sum(weights) / (generated events), with the numbers of saved and generated events; the last event has the final value.

Weights and cross section are proportional to |U|^2 while the branching ratios are small. They use Pythia's meson production: for b hadrons it is about 1.5x above the measured cross section, with too few Lambda_b. The parent of the HNL is in the event record, so the events can be reweighted to data.

Validation against counting the HNLs of events in which the parents decay as in the SM (with their HNL channels): b (m_N = 1 GeV, |U|^2 = 10), c (1 GeV, 1), s (0.3 GeV, 0.1) and pi (0.02 GeV, 0.1) agree to -2.0, +3.6, -3.3 and -1.3 %, within the statistical uncertainties of 2-4 %.

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
