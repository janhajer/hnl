#pragma once

#include "Pythia8/SigmaProcess.h"

namespace hnl {

struct Sigma : public Pythia8::Sigma2Process {
    Sigma(int id_, double coupling_);
protected:
    virtual int code() const override;
    virtual std::string inFlux() const override;
    virtual void initProc() override;
//     virtual int resonanceA() const override;
//     virtual int resonanceB() const override;
    virtual void sigmaKin() override;
    virtual bool convertM2() const override {return true;}
    virtual double sigmaHat() override;
    virtual void setIdColAcol() override;
    virtual double weightDecay(Pythia8::Event& , int , int ) override;
    virtual std::string name() const override;
//     virtual bool isSChannel() const override;

  virtual int id3Mass() const override;
  virtual int id4Mass() const override;
private:
    int id;
    double coupling;
    double mass2;
    double M2;
};

}
