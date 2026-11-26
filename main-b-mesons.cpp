// HNLs from B meson, Lambda_b and tau decays in pp collisions at 14 TeV with Pythia 8.3.
// b-mesons <m_N [GeV]> <|U|^2> <events> <mode> [seed] [flavour: e, mu (default), tau: the neutrino the HNL mixes with]
// Writes b-mesons-<m_N>.csv, one line per HNL from a B meson or tau: weight [pb], parent, HNL four-momentum [GeV], production and decay vertex [mm],
// parent four-momentum [GeV].
// mode 0: SM branching ratios, the weight is the cross section per event.
// mode 1: all B mesons decay into HNLs, each HNL is weighted with the branching ratio of its parent;
//         fine for HNL distributions, but the other B of the event does not decay as in the SM.
// mode 2: one B per event decays into HNLs, picked according to the branching ratios; all other hadrons decay as in the SM.
//         The weight is the sum of the branching ratios of the B mesons in the event.
// mode 3: all hadrons decay as in the SM, every tau decays into HNLs (m_N < m_tau); the weight is BR(tau -> N X).
// Modes 0 and 3 include tau -> N X (tau from D_s, B, ... decays), modes 1 and 2 only B mesons and Lambda_b.
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>

#include "Pythia8/Pythia.h"

#include "BaryonDecays.hh"
#include "MesonDecays.hh"
#include "ResonanceWidths.hh"
#include "TauResonance.hh"
#include "pythia.hh"
#include "string.hh"

using hnl::parse;

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    std::string usage = "usage: b-mesons <m_N [GeV] > 0> <|U|^2 > 0> <events > 0> <mode: 0 SM, 1 all B, 2 one B, 3 all tau> [seed] [flavour: e, mu, tau]";
    double mass = 0., coupling = 0.;
    int events = 0, mode = -1, seed = 0;
    int id_N = arguments.size() < 6 ? hnl::heavy_neutral_muon : arguments.at(5) == "e" ? hnl::heavy_neutral_electron : arguments.at(5) == "mu" ? hnl::heavy_neutral_muon : arguments.at(5) == "tau" ? hnl::heavy_neutral_tau : 0;
    if (arguments.size() < 4 || arguments.size() > 6 || id_N == 0
        || !parse(arguments.at(0), mass) || mass <= 0.
        || !parse(arguments.at(1), coupling) || coupling <= 0.
        || !parse(arguments.at(2), events) || events <= 0
        || !parse(arguments.at(3), mode) || mode < 0 || mode > 3
        || (arguments.size() >= 5 && (!parse(arguments.at(4), seed) || seed < 0 || seed > 900000000))) {
        hnl::print(usage);
        return 1;
    }

    Pythia8::Pythia pythia(hnl::xml_dir(), false);
    hnl::set_pythia_production(pythia);
    hnl::set_pythia_quiet(pythia);
    pythia.readString("SoftQCD:nonDiffractive = on");
    if (arguments.size() > 4) {
        pythia.readString("Random:setSeed = on");
        pythia.readString("Random:seed = " + arguments.at(4));
    }
    std::vector<int> mesons{511, 521, 531, 541, 5122}; // 5122: Lambda_b -> Lambda_c l N
    if (mode == 0 || mode == 3) mesons.push_back(15); // tau -> N X
    pythia.addResonancePtr(std::make_shared<hnl::NeutrinoResonance>(pythia, hnl::neutrino_coupling(coupling), mass, id_N));
    for (auto meson : mesons) {
        if (meson == 15) pythia.addResonancePtr(std::make_shared<hnl::TauResonance>(pythia, hnl::neutrino_coupling(coupling)));
        else if (meson == 5122) pythia.addResonancePtr(std::make_shared<hnl::BaryonResonance>(pythia, hnl::neutrino_coupling(coupling), meson));
        else pythia.addResonancePtr(std::make_shared<hnl::MesonResonance>(pythia, hnl::neutrino_coupling(coupling), meson));
    }
    auto decays = std::make_shared<hnl::MesonDecays>(pythia, mesons);
    decays->force_tau = mode == 3;
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
    if (mode == 3 && branching_ratio[15] <= 0.) {
        hnl::print("no tau -> N X for this mass and flavour");
        return 1;
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
            if (particle.idAbs() != id_N) continue;
            int parent = pythia.event[particle.mother1()].idAbs();
            if (!branching_ratio.count(parent)) continue;
            if (mode == 2 && particle.mother1() != decays->forced) continue;
            if (mode == 3 && parent != 15) continue;
            double weight = mode == 0 ? 1. : mode == 2 ? sum : branching_ratio[parent];
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
