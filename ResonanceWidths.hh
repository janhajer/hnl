#pragma once

#include <functional>
#include "Pythia8/ResonanceWidths.h"
#include "ThreeBodyWidth.hh"

namespace Pythia8{
    class Pythia;
}

namespace hnl
{

// struct ResonanceWidths : public Pythia8::ResonanceWidths {
//     ResonanceWidths(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int id_from);
// protected:
//     virtual bool initBSM() override;
//     virtual bool allowCalc() override;
//     virtual void initConstants() override;
//     virtual void calcPreFac(bool called_from_init = false) override;
//     virtual void calcWidth(bool called_from_init = false) override;
// private:
//     bool can_two_body();
//     bool can_three_body(int id);
//     void add_two_body(Pythia8::ParticleDataEntry& particle);
//     void add_three_body(Pythia8::ParticleDataEntry& particle, int id);
//     bool getChannels();
// private:
//     std::function<double (int id_heavy, int id_light)> neutrino_coupling;
//     ThreeBodyWidth three_body_width;
// };

struct MesonResonance : public Pythia8::ResonanceWidths {
    MesonResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int id_from);
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
    double sum = 0.;
};

struct NeutrinoResonance : public Pythia8::ResonanceWidths {
    NeutrinoResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int id_from);
protected:
    virtual bool initBSM() override;
    virtual bool allowCalc() override;
    virtual void initConstants() override;
    virtual void calcPreFac(bool called_from_init = false) override;
    virtual void calcWidth(bool called_from_init = false) override;
private:
    double get_mass(int id);
    bool can_two_body(int meson);
    bool can_three_body();
    void add_three_body();
    void add_two_body(int id);
    std::vector<int> mesons();
    double CKM2(int meson);
    double CKM2(int up, int down);
    double NW(int up, int down);
    double NZ(int id2Abs);
    double Cf1(int neutrino, int fermion);
    double Cf2(int neutrino, int fermion);
    double correction_factor(int id);
private:
    std::function<double (int id_heavy, int id_light)> neutrino_coupling;
    NeutrinoThreeBodyWidth three_body_width;
    double sum = 0.;
};

}
