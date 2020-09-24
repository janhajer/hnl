#pragma once

#include "Pythia8/Pythia.h"

#include "id.hh"
#include "range.hh"
#include "math.hh"

namespace hnl {

inline void set_pythia_init_quiet(Pythia8::Pythia& pythia) {
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showMultipartonInteractions = off");
}

inline void set_pythia_next_quiet(Pythia8::Pythia& pythia) {
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberCount = 0");
}

inline void set_pythia_quiet(Pythia8::Pythia& pythia) {
    set_pythia_init_quiet(pythia);
    set_pythia_next_quiet(pythia);
}

inline void set_pythia_production(Pythia8::Pythia& pythia) {
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000.");
    pythia.readString("ResonanceWidths:minWidth = 1E-30");
}

inline void set_pythia_passive(Pythia8::Pythia& pythia) {
    set_pythia_init_quiet(pythia);
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("ResonanceWidths:minWidth = 1E-30");
}

inline void set_pythia_branching_ratios(Pythia8::Pythia& pythia) {
    set_pythia_passive(pythia);
    set_pythia_production(pythia);
}

inline void set_pythia_read_hepmc(Pythia8::Pythia& pythia) {
    set_pythia_passive(pythia);
    set_pythia_next_quiet(pythia);
}

inline void set_pythia_read_lhe(Pythia8::Pythia& pythia, std::string const& path) {
    set_pythia_quiet(pythia);
    pythia.readString("Next:numberShowLHA = 0");
    pythia.readString("ResonanceWidths:minWidth = 1E-30");
    pythia.readString("SLHA:readFrom = 0");
    pythia.readString("LesHouches:setLifetime = 2");
    pythia.readString("Beams:frameType = 4"); // LHEF
    pythia.readString("Beams:LHEF = " + path);
}

inline void set_pythia_sigma(Pythia8::Pythia& pythia) {
    set_pythia_production(pythia);
    set_pythia_quiet(pythia);
}

inline void set_pythia_stable(Pythia8::Pythia& pythia, int id, double mass) {
    pythia.particleData.m0(id, mass);
    pythia.particleData.mMin(id, mass / 2);
    pythia.particleData.mMax(id, mass * 2);
    pythia.particleData.particleDataEntryPtr(id)->clearChannels();
}

inline void set_pythia_mesons(Pythia8::Pythia& pythia, int id, double mass) {
    set_pythia_stable(pythia, id, mass);
    set_pythia_production(pythia);
    set_pythia_next_quiet(pythia);
    pythia.readString("Main:numberOfEvents = 100000");
//     pythia.readString("Main:numberOfEvents = 1000");
}

inline void set_pythia_write_hepmc(Pythia8::Pythia& pythia, int id, double mass) {
    set_pythia_mesons(pythia, id, mass);
    pythia.readString("Bottomonium:all = on");
    if (mass < 1.96849) pythia.readString("Charmonium:all = on");
    if (mass < .493677) pythia.readString("HardQCD:all = on");
//     pythia.readString("Onia:all(3S1) = on");
//     pythia.readString("SoftQCD:nonDiffractive = on");
//     pythia.readString("HardQCD:qqbar2ccbar  = on");
//     pythia.readString("HardQCD:gg2bbbar  = on");
//     pythia.readString("PhaseSpace:pTHatMin = .1");
}

inline void set_pythia_minimum_bias(Pythia8::Pythia& pythia, int id, double mass) {
    set_pythia_mesons(pythia, id, mass);
    pythia.readString("SoftQCD:nonDiffractive = on");
}

inline auto for_each(Pythia8::ParticleDataEntry const& particle, std::function<void(Pythia8::DecayChannel const&)> const& function) {
    for (auto channel_number : irange(particle.sizeChannels())) function(particle.channel(channel_number));
}

inline auto for_each(Pythia8::ParticleDataEntry& particle, std::function<void(Pythia8::DecayChannel&)> const& function) {
    for (auto channel_number : irange(particle.sizeChannels())) function(particle.channel(channel_number));
}

inline auto for_each(Pythia8::DecayChannel const& channel, std::function<void(int product)> const& function) {
    for (auto pos : irange(channel.multiplicity())) function(channel.product(pos));
}

inline auto for_each_until(Pythia8::DecayChannel const& channel, std::function<bool(int product)> const& function) {
    for (auto pos : irange(channel.multiplicity())) if(function(channel.product(pos))) return true;
    return false;
}

inline bool has_neutrino(Pythia8::DecayChannel const& channel) {
//     for (auto heavy : heavy_neutral_leptons()) if(for_each_until(channel, [heavy](int product) {
//         if (product == heavy) return true;
//     })) return true;
    for (auto heavy : heavy_neutral_leptons()) for(auto mult : irange(channel.multiplicity())) if(channel.product(mult) == heavy) return true;
    return false;
}

inline double tau_to_Gamma(double tau) { //mm/c->GeV
    return 1.97327E-13 / tau;
}

inline double Gamma_to_tau(double Gamma) { //GeV -> mm/c
    return 1.97327E-13 / Gamma;
}

}

