#pragma once

#include <functional>
#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/StandardModel.h"
#include "ThreeBodyWidth.hh"

namespace Pythia8{
    class Pythia;
}

namespace hnl
{

struct MesonResonance : public Pythia8::ResonanceWidths {
    MesonResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int id_from);
    Pythia8::ParticleDataEntry & calculate_widths();
    bool init(Pythia8::Info* infoPtrIn) override;
protected:
    virtual bool initBSM() override;
    virtual bool allowCalc() override;
    virtual void initConstants() override;
    virtual void calcPreFac(bool called_from_init = false) override;
    virtual void calcWidth(bool called_from_init = false) override;
private:
    bool can_two_body();
    bool can_three_body(int meson);
    void add_two_body();
    void add_three_body(int meson);
    std::vector<int> mesons();
    double CKM2(int id);
    double CKM2(int id_1, int id_2);
private:
    std::function<double (int id_heavy, int id_light)> neutrino_coupling;
    MesonThreeBodyWidth three_body_width;
    Pythia8::ParticleDataEntryPtr particlePtr; // the base class only keeps a weak_ptr
    double sum = 0.;
    bool channels_added = false;
    double tau0_sm;
    double width_sm;
};

struct NeutrinoResonance : public Pythia8::ResonanceWidths {
    // majorana: also the charge conjugated decays, which doubles the width; otherwise only N -> l- X+ and N -> nu X
    NeutrinoResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& coupling, double mass, int id, bool majorana_ = true);
    Pythia8::ParticleDataEntry & calculate_widths();
protected:
    virtual bool initBSM() override;
    virtual bool allowCalc() override;
    virtual void initConstants() override;
    virtual void calcPreFac(bool called_from_init = false) override;
    virtual void calcWidth(bool called_from_init = false) override;
private:
    double get_mass(int id);
    void add_charge_conjugates();
    double CKM2(int up, int down);
    double NZ(int id2Abs);
    double Cf1(int neutrino, int fermion);
    double Cf2(int neutrino, int fermion);
    double correction_factor(int id);
    double coupling(int neutrino);
    double qcd_factor();
    double width_charged_meson(int lepton, int meson);
    double width_neutral_meson(int neutrino, int meson);
    double width_charged_fermions(int lepton, int up, int down);
    double width_neutral_fermions(int neutrino, int fermion);
    double width_invisible(int neutrino, int neutrino_2);
    double remainder_charged(int lepton, int up, int down);
    double remainder_neutral(int neutrino, int quark);
private:
    std::function<double (int id_heavy, int id_light)> neutrino_coupling;
    NeutrinoThreeBodyWidth three_body_width;
    Pythia8::AlphaStrong alpha_strong;
    bool majorana;
    Pythia8::ParticleDataEntryPtr particlePtr; // the base class only keeps a weak_ptr
};

}
