#pragma once

#include <functional>
#include "Pythia8/SigmaProcess.h"

namespace hnl {

struct Sigma : public Pythia8::Sigma2Process {
    Sigma(std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int heavy_, int light_);
protected:
    virtual int code() const override;
    virtual std::string inFlux() const override;
    virtual void initProc() override;
//     virtual int resonanceA() const override {return 23;};
//     virtual int resonanceB() const override;
    virtual void sigmaKin() override;
    virtual bool convertM2() const override {
        return true;
    }
    virtual double sigmaHat() override;
    virtual void setIdColAcol() override;
    virtual double weightDecay(Pythia8::Event&, int, int ) override;
    virtual std::string name() const override;
//     virtual bool isSChannel() const override {return true;};

    virtual int id3Mass() const override;
    virtual int id4Mass() const override;
private:
    int heavy;
    int light;
    std::function<double (int id_heavy, int id_light)> neutrino_coupling;
    double mass2;
    double M2;
};

}
