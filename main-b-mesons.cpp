// HNLs from B meson and Lambda_b decays in pp collisions at 14 TeV with Pythia 8.3.
// b-mesons <m_N [GeV]> <|U_mu|^2> <events> <mode> [seed]
// Writes b-mesons-<m_N>.csv, one line per HNL from a B meson: weight [pb], parent, HNL four-momentum [GeV], production and decay vertex [mm],
// parent four-momentum [GeV].
// mode 0: SM branching ratios, the weight is the cross section per event.
// mode 1: all B mesons decay into HNLs, each HNL is weighted with the branching ratio of its parent;
//         fine for HNL distributions, but the other B of the event does not decay as in the SM.
// mode 2: one B per event decays into HNLs, picked according to the branching ratios; all other hadrons decay as in the SM.
//         The weight is the sum of the branching ratios of the B mesons in the event.
#include <fstream>
#include <iostream>
#include <map>
#include <memory>

#include "Pythia8/Pythia.h"

#include "BaryonDecays.hh"
#include "MesonDecays.hh"
#include "ResonanceWidths.hh"
#include "pythia.hh"
#include "string.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    if (arguments.size() < 4) {
        hnl::print("usage: b-mesons <m_N [GeV]> <|U_mu|^2> <events> <mode: 0 SM, 1 all B, 2 one B> [seed]");
        return 1;
    }
    double mass = hnl::to_double(arguments.at(0));
    double coupling = hnl::to_double(arguments.at(1));
    int events = std::stoi(arguments.at(2));
    int mode = std::stoi(arguments.at(3));

    Pythia8::Pythia pythia(hnl::xml_dir(), false);
    hnl::set_pythia_production(pythia);
    hnl::set_pythia_quiet(pythia);
    pythia.readString("SoftQCD:nonDiffractive = on");
    if (arguments.size() > 4) {
        pythia.readString("Random:setSeed = on");
        pythia.readString("Random:seed = " + arguments.at(4));
    }
    std::vector<int> mesons{511, 521, 531, 541, 5122}; // 5122: Lambda_b -> Lambda_c l N
    pythia.addResonancePtr(std::make_shared<hnl::NeutrinoResonance>(pythia, hnl::neutrino_coupling(coupling), mass, hnl::heavy_neutral_muon));
    for (auto meson : mesons) {
        if (meson == 5122) pythia.addResonancePtr(std::make_shared<hnl::BaryonResonance>(pythia, hnl::neutrino_coupling(coupling), meson));
        else pythia.addResonancePtr(std::make_shared<hnl::MesonResonance>(pythia, hnl::neutrino_coupling(coupling), meson));
    }
    auto decays = std::make_shared<hnl::MesonDecays>(pythia, mesons);
    pythia.setDecayPtr(decays);
    if (!pythia.init()) return 1;

    std::map<int, double> branching_ratio; // into HNLs
    for (auto meson : mesons) {
        auto& entry = *pythia.particleData.particleDataEntryPtr(meson);
        hnl::for_each(entry, [&](Pythia8::DecayChannel& channel) {
            if (hnl::has_neutrino(channel)) branching_ratio[meson] += channel.bRatio();
            else if (mode == 1) channel.onMode(0);
        });
        if (mode == 2) entry.setMayDecay(false); // decayed after the choice below
        hnl::print("BR(", meson, "-> N X) =", branching_ratio[meson]);
    }

    std::ofstream file("b-mesons-" + arguments.at(0) + ".csv");
    file << "weight_pb,parent,E,px,py,pz,x,y,z,t,x_dec,y_dec,z_dec,t_dec,parent_E,parent_px,parent_py,parent_pz\n";
    for (int event = 0; event < events; ++event) {
        if (!pythia.next()) continue;
        double sum = 0.;
        if (mode == 2) {
            std::vector<std::pair<int, double>> candidates;
            for (int line = 0; line < pythia.event.size(); ++line) if (pythia.event[line].isFinal() && branching_ratio.count(pythia.event[line].idAbs())) {
                sum += branching_ratio[pythia.event[line].idAbs()];
                candidates.emplace_back(line, sum);
            }
            if (candidates.empty()) continue;
            double pick = sum * pythia.rndm.flat();
            for (auto [line, cumulative] : candidates) if (pick < cumulative) {
                decays->forced = line;
                break;
            }
            for (auto meson : mesons) pythia.particleData.particleDataEntryPtr(meson)->setMayDecay(true);
            pythia.forceHadronLevel(false); // moreDecays() can hang on HNL decays into quarks
            for (auto meson : mesons) pythia.particleData.particleDataEntryPtr(meson)->setMayDecay(false);
        }
        for (int line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            if (particle.idAbs() != hnl::heavy_neutral_muon) continue;
            int parent = pythia.event[particle.mother1()].idAbs();
            if (!branching_ratio.count(parent)) continue;
            if (mode == 2 && particle.mother1() != decays->forced) continue;
            double weight = mode == 0 ? 1. : mode == 1 ? branching_ratio[parent] : sum;
            auto const& p = particle.p();
            auto const& v = particle.vProd();
            auto const& d = particle.vDec();
            auto const& q = pythia.event[particle.mother1()].p();
            file << weight << ',' << parent << ',' << p.e() << ',' << p.px() << ',' << p.py() << ',' << p.pz() << ',' << v.px() << ',' << v.py() << ',' << v.pz() << ',' << v.e() << ',' << d.px() << ',' << d.py() << ',' << d.pz() << ',' << d.e() << ',' << q.e() << ',' << q.px() << ',' << q.py() << ',' << q.pz() << '\n';
        }
        decays->forced = 0;
    }
    // the weights above are per event; multiply by sigma / events
    file << "#sigma_pb_per_event," << pythia.info.sigmaGen() * 1e9 / pythia.info.nAccepted() << '\n';
}
