#pragma once

#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/SusyCouplings.h"
#include "Pythia8/Settings.h"
#include "Pythia8/Pythia.h"
#include "ThreeBodyWidth.hh"

namespace neutrino
{


struct BRatio {
    int onMode;
    double bRatio;
    int meMode;
};

struct MesonResonance : public Pythia8::ResonanceWidths {
    MesonResonance(Pythia8::Pythia & pythia, double neutrino_coupling_, int id_from_);
    void AddMissingChannels(Pythia8::ParticleData & particle_data);
protected:
    virtual bool initBSM() override;
    virtual bool allowCalc() override;
    virtual void initConstants() override;
    virtual void calcPreFac(bool = false) override;
    virtual void calcWidth(bool calledFromInit = false) override;
private:
    bool can_two_body();
    bool can_three_body(int id);
    void add_two_body(Pythia8::ParticleDataEntry & particle, int neutrino);
    void add_three_body(Pythia8::ParticleDataEntry & particle, int neutrino, int id);
    bool getChannels();
    double CKM2(int id);
    double CKM2(int id_1,int id_2);
    double NeutU(int id_heavy, int id_light);
    double neutrino_coupling;
    double width;
    ThreeBodyWidth three_body_width;
    Pythia8::CoupSM standard_model;
    Pythia8::ParticleDataEntry particle_data_entry;
    std::map<std::tuple<int, int, int, int, int>, BRatio> channels;
};

}
