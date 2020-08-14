#pragma once

#include "Pythia8/Pythia.h"

#include "math.hh"
#include "id.hh"
#include "io.hh"

namespace hnl {

namespace {

const bool debug = false;

}

void set_pythia_production(Pythia8::Pythia& pythia) {
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000.");
    pythia.readString("ResonanceWidths:minWidth = 1E-30");
}

void set_pythia_init(Pythia8::Pythia& pythia) {
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showMultipartonInteractions = off");
}

void set_pythia_passive(Pythia8::Pythia& pythia) {
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("ResonanceWidths:minWidth = 1E-30");
}

void set_pythia_next(Pythia8::Pythia& pythia) {
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowInfo = 0");
}


void set_pythia_read_hepmc(Pythia8::Pythia& pythia) {
    set_pythia_init(pythia);
    set_pythia_next(pythia);
    set_pythia_passive(pythia);
}

void set_pythia_branching_ratios(Pythia8::Pythia& pythia) {
    set_pythia_production(pythia);
    set_pythia_init(pythia);
    set_pythia_passive(pythia);
}

void set_pythia_sigma(Pythia8::Pythia& pythia) {
    set_pythia_production(pythia);
    set_pythia_init(pythia);
    set_pythia_next(pythia);
}

auto has_neutrino = [](Pythia8::DecayChannel const& channel) {
    for (auto heavy : heavy_neutral_leptons()) {
        if (channel.product(0) == heavy || channel.product(1) == heavy || channel.product(2) == heavy || channel.product(3) == heavy || channel.product(4) == heavy) return true;
    }
    return false;
};

struct Loop {
    Loop(double min, int steps_) : m_min(min), steps(steps_) {}
    double m_min;
    int steps;
    double mass(double max, int step) const {
        return log_scale(m_min, max, step, steps);
    }
};

void set_pythia_write_hepmc(Pythia8::Pythia& pythia, double mass) {
    pythia.particleData.m0(heavy_neutrino, mass);
    pythia.particleData.mMin(heavy_neutrino, 0.);
    pythia.particleData.particleDataEntryPtr(heavy_neutrino)->clearChannels();
    set_pythia_production(pythia);
    set_pythia_next(pythia);
    pythia.readString("Bottomonium:all = on");
    if (mass < 1.96849) pythia.readString("Charmonium:all = on");
    if (mass < .493677) pythia.readString("HardQCD:all = on");
//     pythia.readString("Onia:all(3S1) = on");
//     pythia.readString("SoftQCD:nonDiffractive = on");
//     pythia.readString("HardQCD:qqbar2ccbar  = on");
//     pythia.readString("HardQCD:gg2bbbar  = on");
//     pythia.readString("PhaseSpace:pTHatMin = .1");
    pythia.readString("Main:numberOfEvents = 100000");
}

}

