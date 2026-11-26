#pragma once

// HNLs from meson decays in inclusive pp collisions (SoftQCD:nonDiffractive) at 14 TeV, written as HepMC3 ascii.
//
// The HNL is nu_Rmu (9900014) with mass m_N, mixing |U_mu|^2 with nu_mu, and stable. The parents are the hadrons of one
// flavour (parents()). Each event is generated with the SM, then one of its parents decays into an HNL: it is picked
// according to the branching ratios into HNLs of all parents of the event and decays with the matrix element of
// MesonDecays; all other hadrons decay as in the SM. A parent is only found if no parent decays into it, so the flavours
// are generated separately: b hadrons decay normally in the c run, and their D mesons are parents there. The event
// weight is the sum of these branching ratios, the probability that the event contains an HNL (for small branching
// ratios). Saved are events whose HNL has 1.4 < |eta| < 3.5.
//
// Each saved event carries the cross section [pb] of HNLs in the acceptance so far, sigma_gen * sum(weights) / (generated
// events), with the numbers of saved and generated events; the last event has the final value. The weights and the cross
// section are proportional to |U_mu|^2 as long as the branching ratios are small, and so can be rescaled.

#include <cmath>
#include <map>
#include <memory>
#include <string>

#include "BaryonDecays.hh"
#include "MesonDecays.hh"
#include "ResonanceWidths.hh"
#include "hepmc.hh"
#include "string.hh"
#include "pythia.hh"

namespace hnl {

namespace hepmc {

namespace {

constexpr bool debug = false;

}

// the HNL: the other two keep Pythia's mass of 500 GeV and are not produced
constexpr int hnl_id = heavy_neutral_muon;

// the parents of one flavour: b, c, s or pi; none of them comes from the decay of another (except B_c -> B_s, rare)
inline std::vector<int> parents(std::string const& flavour) {
    if (flavour == "b") return {511, 521, 531, 541, 5122, 553};
    if (flavour == "c") return {411, 421, 431, 443};
    if (flavour == "s") return {321, 130, 310};
    if (flavour == "pi") return {211};
    return {};
}

// HNL in the event that comes from the meson in line `forced`, 0 if none
inline int hnl_from(Pythia8::Event const& event, int forced) {
    for (int line = 0; line < event.size(); ++line) if (event[line].idAbs() == hnl_id && event[line].mother1() == forced) return line;
    return 0;
}

inline void write(Pythia8::Pythia& pythia, MesonDecays& decays, std::map<int, double> const& branching_ratio, double mass, double coupling, std::string const& flavour, int events) {
    auto run_info = std::make_shared<HepMC3::GenRunInfo>();
    run_info->tools().push_back({"Pythia8", std::to_string(pythia.parm("Pythia:versionNumber")), "pp at " + std::to_string(pythia.info.eCM()) + " GeV"});
    run_info->set_weight_names(pythia.info.weightNameVector());
    run_info->add_attribute("hnl_id", std::make_shared<HepMC3::IntAttribute>(hnl_id));
    run_info->add_attribute("hnl_mass_GeV", std::make_shared<HepMC3::DoubleAttribute>(mass));
    run_info->add_attribute("hnl_U2", std::make_shared<HepMC3::DoubleAttribute>(coupling));
    run_info->add_attribute("parent_flavour", std::make_shared<HepMC3::StringAttribute>(flavour));
    run_info->add_attribute("weight", std::make_shared<HepMC3::StringAttribute>("sum of the branching ratios into HNLs of the mesons of the event"));
    HepMC3::WriterAscii file(std::to_string(mass) + "-" + flavour + ".hepmc", run_info);
    HepMC3::Pythia8ToHepMC3 converter;

    long generated = 0; // events with pythia.next() == true
    int saved = 0;
    double sum = 0.; // of the weights of the saved events
    double sum2 = 0.;
    auto cross_section = [&]() -> std::pair<double, double> { // pb
        double sigma = pythia.info.sigmaGen() * 1e9 * sum / generated;
        double error = sum > 0. ? sigma * std::sqrt(sqr(pythia.info.sigmaErr() / pythia.info.sigmaGen()) + sum2 / sqr(sum)) : 0.;
        return {sigma, error};
    };

    while (saved < events) {
        if (!pythia.next()) continue;
        ++generated;
        double weight = 0.;
        std::vector<std::pair<int, double>> candidates;
        for (int line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            if (!particle.isFinal() || !branching_ratio.count(particle.idAbs())) continue;
            weight += branching_ratio.at(particle.idAbs());
            candidates.emplace_back(line, weight);
        }
        if (weight <= 0.) continue;
        double pick = weight * pythia.rndm.flat();
        for (auto [line, cumulative] : candidates) if (pick < cumulative) {
            decays.forced = line;
            break;
        }
        for (auto const& [meson, ratio] : branching_ratio) pythia.particleData.particleDataEntryPtr(meson)->setMayDecay(true);
        bool decayed = pythia.forceHadronLevel(false);
        for (auto const& [meson, ratio] : branching_ratio) pythia.particleData.particleDataEntryPtr(meson)->setMayDecay(false);
        int line = decayed ? hnl_from(pythia.event, decays.forced) : 0;
        decays.forced = 0;
        if (line == 0) {
            if (debug) print("no HNL after the forced decay in event", generated);
            continue;
        }
        double eta = std::abs(pythia.event[line].eta());
        if (!(1.4 < eta && eta < 3.5)) continue;

        ++saved;
        sum += weight;
        sum2 += sqr(weight);
        if (debug) print("saved event", saved, "of", generated, "HNL from", pythia.event[pythia.event[line].mother1()].id(), "eta", eta, "weight", weight);
        HepMC3::GenEvent event(run_info, HepMC3::Units::GEV, HepMC3::Units::MM);
        converter.fill_next_event(pythia, event);
        for (auto& value : event.weights()) value *= weight;
        auto [sigma, error] = cross_section();
        auto cross_section_ptr = std::make_shared<HepMC3::GenCrossSection>();
        event.set_cross_section(cross_section_ptr);
        cross_section_ptr->set_cross_section(sigma, error, saved, generated);
        file.write_event(event);
    }
    file.close();

    pythia.stat();
    auto [sigma, error] = cross_section();
    print("Saved", saved, "of", generated, "events; cross section of HNLs in the acceptance", sigma, "+-", error, "pb at |U|^2 =", coupling);
}

// write-hepmc <m_N [GeV]> <|U_mu|^2> <b|c|s|pi> [events] [seed]
inline bool write_hepmc(double mass, double coupling, std::string const& flavour, int events = 1000, int seed = -1) {
    print("Generating events for HNLs with mass", mass, "GeV and |U|^2 =", coupling);
    auto mesons = parents(flavour);
    if (mesons.empty()) {
        print("Unknown parent flavour", flavour);
        return false;
    }
    Pythia8::Pythia pythia(xml_dir(), false);
    set_pythia_minimum_bias(pythia, hnl_id, mass);
    if (seed >= 0) {
        pythia.readString("Random:setSeed = on");
        pythia.readString("Random:seed = " + std::to_string(seed));
    }
    for (auto meson : mesons) {
        if (meson == 5122) pythia.addResonancePtr(std::make_shared<BaryonResonance>(pythia, neutrino_coupling(coupling), meson));
        else pythia.addResonancePtr(std::make_shared<MesonResonance>(pythia, neutrino_coupling(coupling), meson));
    }
    auto decays = std::make_shared<MesonDecays>(pythia, mesons);
    pythia.setDecayPtr(decays);
    if (!pythia.init()) return false;

    std::map<int, double> branching_ratio; // into HNLs, of the mesons that have such decays
    for (auto meson : mesons) {
        double ratio = 0.;
        for_each(*pythia.particleData.particleDataEntryPtr(meson), [&ratio](Pythia8::DecayChannel const& channel) {
            if (has_neutrino(channel)) ratio += channel.bRatio();
        });
        print("BR(", meson, "-> N X) =", ratio);
        if (ratio <= 0.) continue;
        branching_ratio[meson] = ratio;
        pythia.particleData.particleDataEntryPtr(meson)->setMayDecay(false); // decayed after the choice in write()
    }
    if (branching_ratio.empty()) {
        print("No meson decays into HNLs of mass", mass, "GeV");
        return false;
    }
    write(pythia, *decays, branching_ratio, mass, coupling, flavour, events);
    return true;
}

}

}
