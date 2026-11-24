#pragma once

#include <map>
#include <tuple>

#include "Pythia8/ParticleDecays.h"
#include "ThreeBodyWidth.hh"

namespace Pythia8{
    class Pythia;
}

namespace hnl
{

// Meson (and Lambda_b) decays into HNLs with the matrix element instead of phase space.
// The HNL channels stay in the decay tables: the channel is picked here as Pythia would, SM channels go back to Pythia.
struct MesonDecays : public Pythia8::DecayHandler {
    MesonDecays(Pythia8::Pythia& pythia, std::vector<int> const& mesons_);
    bool decay(std::vector<int>& ids, std::vector<double>& masses, std::vector<Pythia8::Vec4>& momenta, int line, Pythia8::Event const&) override;
    std::vector<int> handledParticles() override;
    // event line of a meson that has to decay into HNLs (0: none)
    int forced = 0;
    // |M|^2 / (4 G_F^2 |V|^2 |U|^2) of meson(p) -> daughter(p_to) lepton(k_l) N(k_N)
    double matrix_element(int id_from, int id_to, Pythia8::Vec4 const& p, Pythia8::Vec4 const& p_to, Pythia8::Vec4 const& k_l, Pythia8::Vec4 const& k_N);
    // density in (q2, cos_theta, phi) at rest; sets the momenta
    double weight(int id_from, int id_to, double m_from, double m_to, double m_l, double m_N, double q2, double cos_theta, double phi, Pythia8::Vec4& p_to, Pythia8::Vec4& k_l, Pythia8::Vec4& k_N);
    // meson at rest -> daughter, lepton and N, distributed according to the matrix element
    bool three_body(int id_from, int id_to, double m_from, double m_to, double m_l, double m_N, Pythia8::Vec4& p_to, Pythia8::Vec4& k_l, Pythia8::Vec4& k_N);
private:
    void restore();
    Pythia8::ParticleData* particle_data;
    Pythia8::Rndm* rndm;
    std::vector<int> mesons;
    MesonThreeBodyWidth three_body_width;
    std::map<std::tuple<int, int, int, int, int>, double> max_weight;
    std::vector<std::tuple<Pythia8::ParticleDataEntryPtr, int, int>> switched_off;
};

}
