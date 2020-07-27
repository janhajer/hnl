#pragma once

#include <functional>
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

// struct Resonance : public Pythia8::ResonanceWidths {
//     Resonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int id_from);
//     void AddMissingChannels(Pythia8::ParticleData& particle_data);
// protected:
//     virtual bool initBSM() override;
//     virtual bool allowCalc() override;
//     virtual void initConstants() override;
//     virtual void calcPreFac(bool calledFromInit = false) override;
//     virtual void calcWidth(bool calledFromInit = false) override;
// private:
//     bool can_two_body();
//     bool can_three_body(int id);
//     void add_two_body(Pythia8::ParticleDataEntry& particle);
//     void add_three_body(Pythia8::ParticleDataEntry& particle, int id);
//     bool getChannels();
//     double CKM2(int id);
//     double CKM2(int id_1, int id_2);
// private:
//     std::function<double (int id_heavy, int id_light)> neutrino_coupling;
//     ThreeBodyWidth three_body_width;
//     Pythia8::CoupSM standard_model;
//     Pythia8::ParticleDataEntry particle_data_entry;
// };

struct MesonResonance : public Pythia8::ResonanceWidths {
    MesonResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int id_from);
protected:
    virtual bool initBSM() override;
    virtual bool allowCalc() override;
    virtual void initConstants() override;
    virtual void calcPreFac(bool calledFromInit = false) override;
    virtual void calcWidth(bool calledFromInit = false) override;
private:
    bool can_two_body();
    bool can_three_body(int id);
    void add_two_body();
    void add_three_body(int id);
    std::vector<int> mesons();
    bool getChannels();
    double CKM2(int id);
    double CKM2(int id_1, int id_2);
private:
    std::function<double (int id_heavy, int id_light)> neutrino_coupling;
    ThreeBodyWidth three_body_width;
};

struct NeutrinoResonance : public Pythia8::ResonanceWidths {
    NeutrinoResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int id_from);
protected:
    virtual bool initBSM() override;
    virtual bool allowCalc() override;
    virtual void initConstants() override;
    virtual void calcPreFac(bool calledFromInit = false) override;
    virtual void calcWidth(bool calledFromInit = false) override;
private:
    bool can_two_body(int id);
    bool can_three_body();
    void add_three_body();
    void add_two_body(int id);
    std::vector<int> mesons();
    bool getChannels();
    double CKM2(int id);
    double CKM2(int up, int down);
    double NW(int up, int down);
    double NZ(int id2Abs);
    double Cf1(int id2Abs,int id3Abs);
    double Cf2(int id2Abs,int id3Abs);
    double correction_factor(int id);
private:
    std::function<double (int id_heavy, int id_light)> neutrino_coupling;
    NeutrinoThreeBodyWidth three_body_width;
    double sum = 0;
};

}
