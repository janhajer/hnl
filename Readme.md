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

`b-mesons <m_N> <|U_mu|^2> <events> <forced>` (`main-b-mesons.cpp`) shows the setup:

* `NeutrinoResonance`: HNL mass, widths and decays (Majorana by default),
* `MesonResonance` for each meson: HNL channels added to its decay table, widths from arXiv:1805.08567,
* `MesonDecays`: decay handler that generates the HNL channels with the matrix element instead of phase space.
