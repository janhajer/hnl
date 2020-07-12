#pragma once

#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/SusyCouplings.h"
#include "Pythia8/Settings.h"
#include "ThreeBodyWidth.hh"

namespace neutrino
{

struct MesonResonance : public Pythia8::ResonanceWidths {
    MesonResonance(Pythia8::Settings& settings, Pythia8::Rndm* rndmPtrIn, double neutrino_coupling_, int idResIn) :
        neutrino_coupling(neutrino_coupling_),
        pseudo_scalar()
    {
        initBasic(idResIn);
       standard_model.init(settings, rndmPtrIn);
    }
private:
    virtual bool allowCalc() override;
    virtual bool initBSM() override;
    bool getChannels(int idPDG); // Locally stored properties and couplings.
    virtual void initConstants() override; // Initialize constants.
    virtual void calcPreFac(bool = false) override; // Calculate various common prefactors for the current mass.
    virtual void calcWidth(bool calledFromInit = false) override; // Calculate width for currently considered channel.
    double CKM2(int id);
    double NeutU(int id_heavy, int id_light);
    double neutrino_coupling;
//     double gf2;
    int id;
    // Three-body stau decaywidth classes
    ThreeBodyWidth pseudo_scalar;
    Pythia8::CoupSM standard_model;
};

}
