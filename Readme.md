# Heavy Neutral Leptons in Pythia8

## Dependencies

* Pythia 8.3 (tested with 8.318)
* Boost headers (Range, Optional, Algorithm), e.g. `libboost-dev`
* CMake >= 3.15
* HepMC3, optional, only for `write-hepmc`, e.g. `libhepmc3-dev`
* Catch2 v3, only for the tests; downloaded at configure time if not installed (e.g. `catch2`)

## Build

```
mkdir build && cd build
cmake .. -DPYTHIA8_DIR=<pythia prefix>   # default: $PYTHIA8, else ~/Development/pythia
make -j
```

`write-hepmc` is built if CMake finds HepMC3; otherwise add `-DHepMC3_DIR=<prefix>/share/HepMC3/cmake`.

## Tests

`ctest` in `build/` runs `tests/`: the form factors against values published with their fits, the kinematic constraints at q^2 = 0, and the SM limit of the widths (HNL of negligible mass, |U|^2 = 1) against measured branching ratios and lattice QCD rates, and the B -> D(*) l N widths of a heavy HNL against their zero-recoil limit (unpublished calculation by S. Biondini and J. Hajer). `-DBUILD_TESTING=OFF` skips them.

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

`b-mesons <m_N> <|U|^2> <events> <mode> [seed] [e|mu|tau]` (`main-b-mesons.cpp`; mode 0: SM branching ratios, 1: all B into HNLs, 2: one B per event into HNLs, 3: every tau into HNLs with all hadrons decaying as in the SM; the last argument is the neutrino the HNL mixes with, mu by default) shows the setup:

* `NeutrinoResonance`: HNL mass, widths and decays (Majorana by default),
* `MesonResonance` for each meson: HNL channels added to its decay table, widths from arXiv:1805.08567 with the inputs below,
* `MesonDecays`: decay handler that generates the HNL channels with the matrix element instead of phase space.
* `BaryonResonance`: Lambda_b -> Lambda_c l N with the lattice QCD form factors of arXiv:1503.01421; |M|^2 from Dirac traces, for any lepton and HNL masses. In the SM limit it gives their Gamma / |V_cb|^2 for e, mu and tau (21.5, 21.5, 7.15 ps^-1).
* `TauResonance`: tau -> N pi, N K, N rho, N K*, N l nubar (mixing with nu_tau) and tau -> nu_tau l N (mixing with nu_e, nu_mu), l = e, mu, widths of arXiv:1805.08567; `MesonDecays` generates the three-body decays with the V-A matrix element. The tau decays as usual otherwise (Pythia's `TauDecays`). Below m_tau this adds the HNLs from taus of D_s and B decays, the main source for mixing with nu_tau.

## HNL decays

Widths of arXiv:1805.08567. Hadronic decays:

* one meson at all masses: `l P` (pi, K, D, D_s, B, B_c), `l V` (rho, K\*, D\*, D_s\*, a_1), `nu P` (pi0, eta, eta', eta_c), `nu V` (rho0, omega, phi, J/psi, a_1). The K\* and a_1 couplings come from tau decays, as g_rho in the paper. The neutral-current factors of omega, phi and J/psi follow from the quark couplings with the paper's decay constants and are smaller than in its table (see `correction_factor`);
* for each quark current, the remainder of the inclusive width (quarks with the QCD correction, alpha_s at three loops) over its one-meson channels, if positive and above the lightest states those channels do not cover. The quark pair fragments in Pythia's hadron level (meMode 94 with the V-A matrix element for the charged current, 93 for the neutral one).

At m_N = m_tau the charged-current channels reproduce tau decays (pi : rho : a_1 : rest).

Decay HNLs with `next()` or `forceHadronLevel()`. `moreDecays()` with its default `allowPartons = false` replaces the quark pair by a random-flavour pair: the charge is not conserved, and below 1.45 GeV `ParticleDecays::pickHadrons` loops forever.

`patches/pythia8318-ParticleDecays.patch` fixes this in Pythia 8.318 (`patch -p1 < .../pythia8318-ParticleDecays.patch` in the Pythia directory, then `make`): `moreDecays(false)` keeps the quark flavours and the flavour loop has a try limit; `next()` is not affected. HNL decays at 0.5-4.5 GeV then conserve the charge and never hang; so do SM decays, where unpatched `moreDecays()` gives the wrong charge in 43% of B+ and 23% of Lambda_b decays. Up to m_N = 2 GeV the multiplicity agrees with `forceHadronLevel()` to 2%; above, `moreDecays()` picks the hadrons by phase space (charged multiplicity +22% at 3 GeV), so `forceHadronLevel()` stays the better choice there.

## Physics inputs

Last checked against the literature: November 2026. Each value is cited again at its place in the code. Inputs that replace Pythia 8.318 defaults are set in `set_physics_inputs()` (`PhysicsInputs.hh`), called by the resonance constructors.

Semileptonic form factors (`FormFactors.cpp`); FLAG 2024 is arXiv:2411.04268, PDG 2025 the 2025 update of the Review of Particle Physics:

| transition | source | notes |
|---|---|---|
| K -> pi | f_+(0) = 0.9698: FLAG 2024 eq. (72); slopes: FlaviaNet, arXiv:1005.2323, table 10 | quadratic f_+, linear f_0; K+ -> pi0 with delta_SU(2) = 0.029 |
| D -> K | FLAG 2024 table 30 (ETM 17D, HPQCD 21A, FNAL/MILC 22) | lattice only |
| D -> pi | FNAL/MILC 22, arXiv:2212.12648, table XIII | FLAG 2024 gives no shape (ETM/FNAL tension) |
| D_s -> eta | LCSR shape of 1805.08567 | normalised to B(D_s -> eta e nu) = 2.27 % (PDG 2025) |
| D -> K* | Melikhov, Stech, hep-ph/0001113 | normalised to B(D+ -> K*0 e nu) = 5.40 % (PDG 2025); no lattice result |
| B -> pi | FLAG 2024 table 50 | lattice and B-factory data; \|V_ub\| of the fit 3.61e-3 |
| B_s -> K | FLAG 2024 table 40 | lattice only |
| B -> D | FLAG 2024 table 52 | lattice and data; \|V_cb\| of the fit 40.0e-3 |
| B_s -> D_s | HPQCD 19, arXiv:1906.00701, table VIII | as FLAG 2024 table 46 |
| B -> D* | FLAG 2024 table 53 (N_f = 2+1) in the BGL form of FNAL/MILC 21, arXiv:2105.14019 | lattice and data; \|V_cb\| of the fit 39.23e-3 |
| B_s -> D_s* | HPQCD 21B, arXiv:2105.11433, table XIII | HPQCD 23 publishes no coefficients |
| B -> rho | BSZ, arXiv:1503.05534, table 14 | LCSR |
| B_s -> K* | BSZ, table 15 | LCSR and lattice |
| B_c -> eta_c | Ivanov, Korner, Santorelli, hep-ph/0007169 | quark model; no lattice result |
| B_c -> J/psi | HPQCD 20, arXiv:2007.06957, table XVI | HPQCD 25 (arXiv:2503.15090) agrees in the total rate |
| Lambda_b -> Lambda_c | Detmold, Lehner, Meinel, arXiv:1503.01421 | still the latest lattice result |

For vector daughters A_0 is shifted by a constant so that A_0(0) = A_3(0) holds exactly (by at most 2.5 %, for B_s -> K*), see `FormFactors.hh`.

Decay constants (`decay_constant()` in `ResonanceWidths.cpp`, conventions of 1805.08567):

| meson | value | source |
|---|---|---|
| pi+, pi0 | 130.2 MeV | FLAG 2024 eq. (79) |
| K+ | 155.7 MeV | FLAG 2024 eq. (80) |
| D+ | 212.0 MeV | FLAG 2024 eq. (129) |
| D_s | 249.9 MeV | FLAG 2024 eq. (130) |
| B+ | 189.4 MeV | FNAL/MILC 17 in FLAG 2024 table 33 |
| B_c | 427 MeV | HPQCD, arXiv:1207.0994 |
| eta, eta' | 81.7, -94.7 MeV | 1805.08567 eqs. (C.11), (C.12) |
| eta_c | 398.1 MeV / sqrt(2) | HPQCD, arXiv:2005.01845 |
| rho | 0.162 GeV^2 | 1805.08567, from tau -> rho nu |
| D*, D_s* | f = 223.5, 268.8 MeV, times the mass | ETMC, arXiv:1707.04529 |
| omega, phi, J/psi | 0.153, 0.234, 1.29 GeV^2 | 1805.08567 table 9 |
| K*, a_1 | 0.1795, 0.2839 GeV^2 | tau branching ratios of PDG 2025 |

Pythia 8.318 defaults replaced in `set_physics_inputs()`:

| quantity | value | source |
|---|---|---|
| \|V_CKM\| | global fit, e.g. V_ub = 3.732e-3, V_cb = 41.83e-3 | PDG 2025 CKM review eq. (12.27) |
| c tau of pi+, K+, K_L, K_S, D, D_s, B, B_s, B_c, Lambda_b | e.g. B_c 0.15289 mm | PDG 2025 |
| branching ratios into tau of D_s, D+, B+, B0 | D_s -> tau nu 5.39 %, B+ (B0) -> tau nu X 2.52 (2.33) %, of which D tau nu 0.86 (0.80) %, D* tau nu 1.47 (1.37) % | PDG 2025 |

Kept from Pythia: masses, V -> e+ e- branching ratios (within 2 % of PDG 2025), sin^2 theta_W. alpha_s(m_Z) = 0.118 in `ResonanceWidths.cpp`.

The global-fit \|V_cb\| and \|V_ub\| are above the exclusive values of the B -> D(*) and B -> pi fits (inclusive-exclusive tension), so these channels are 7-14 % above the measured branching ratios; with the \|V\| of their fits they agree (`tests/widths.cpp`). Pythia limits V_cb to [0.0418, 0.0426].
