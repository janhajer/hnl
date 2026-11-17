// HNLs from B meson decays in pp collisions at 14 TeV with Pythia 8.3.
// b-mesons <m_N [GeV]> <|U_mu|^2> <events> <forced: 0|1> [seed]
// Writes b-mesons-<m_N>.csv, one line per HNL from a B meson: weight [pb], parent, HNL four-momentum [GeV], production and decay vertex [mm].
// forced = 1: all B mesons decay into HNLs and each HNL is weighted with the branching ratio of its parent;
// fine for HNL distributions, but the other B of the event does not decay as in the SM.
// forced = 0: SM branching ratios, the weight is the cross section per event.
#include <fstream>
#include <iostream>
#include <map>
#include <memory>

#include "Pythia8/Pythia.h"

#include "MesonDecays.hh"
#include "ResonanceWidths.hh"
#include "pythia.hh"
#include "string.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    if (arguments.size() < 4) {
        hnl::print("usage: b-mesons <m_N [GeV]> <|U_mu|^2> <events> <forced: 0|1> [seed]");
        return 1;
    }
    double mass = hnl::to_double(arguments.at(0));
    double coupling = hnl::to_double(arguments.at(1));
    int events = std::stoi(arguments.at(2));
    bool forced = std::stoi(arguments.at(3));

    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    hnl::set_pythia_production(pythia);
    hnl::set_pythia_quiet(pythia);
    pythia.readString("SoftQCD:nonDiffractive = on");
    if (arguments.size() > 4) {
        pythia.readString("Random:setSeed = on");
        pythia.readString("Random:seed = " + arguments.at(4));
    }
    std::vector<int> mesons{511, 521, 531, 541};
    pythia.addResonancePtr(std::make_shared<hnl::NeutrinoResonance>(pythia, hnl::neutrino_coupling(coupling), mass, hnl::heavy_neutral_muon));
    for (auto meson : mesons) pythia.addResonancePtr(std::make_shared<hnl::MesonResonance>(pythia, hnl::neutrino_coupling(coupling), meson));
    pythia.setDecayPtr(std::make_shared<hnl::MesonDecays>(pythia, mesons));
    if (!pythia.init()) return 1;

    std::map<int, double> branching_ratio; // into HNLs
    for (auto meson : mesons) {
        auto& entry = *pythia.particleData.particleDataEntryPtr(meson);
        hnl::for_each(entry, [&](Pythia8::DecayChannel& channel) {
            if (hnl::has_neutrino(channel)) branching_ratio[meson] += channel.bRatio();
            else if (forced) channel.onMode(0);
        });
        hnl::print("BR(", meson, "-> N X) =", branching_ratio[meson]);
    }

    std::ofstream file("b-mesons-" + arguments.at(0) + ".csv");
    file << "weight_pb,parent,E,px,py,pz,x,y,z,t,x_dec,y_dec,z_dec,t_dec\n";
    for (int event = 0; event < events; ++event) {
        if (!pythia.next()) continue;
        for (int line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            if (particle.idAbs() != hnl::heavy_neutral_muon) continue;
            int parent = pythia.event[particle.mother1()].idAbs();
            if (!branching_ratio.count(parent)) continue;
            double weight = forced ? branching_ratio[parent] : 1.;
            auto const& p = particle.p();
            auto const& v = particle.vProd();
            auto const& d = particle.vDec();
            file << weight << ',' << parent << ',' << p.e() << ',' << p.px() << ',' << p.py() << ',' << p.pz() << ',' << v.px() << ',' << v.py() << ',' << v.pz() << ',' << v.e() << ',' << d.px() << ',' << d.py() << ',' << d.pz() << ',' << d.e() << '\n';
        }
    }
    // the weights above are per event; multiply by sigma / events
    file << "#sigma_pb_per_event," << pythia.info.sigmaGen() * 1e9 / pythia.info.nAccepted() << '\n';
}
