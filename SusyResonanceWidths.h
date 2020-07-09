#pragma once

#include "Pythia8/ResonanceWidths.h"
#include "Pythia8/SusyCouplings.h"
#include "SusyWidthFunctions.h"

namespace neutrino {

class SUSYResonanceWidths : public Pythia8::ResonanceWidths{
public:
  SUSYResonanceWidths() {}
protected:
  virtual bool initBSM();
  virtual bool allowCalc();
  virtual bool getChannels(int) { return false; };

  double integrateGauss( WidthFunction* widthFn, double, double, double);

  Pythia8::CoupSM* standard_model;
};

// The ResonanceSlepton class handles the Slepton/Sneutrino resonances.

class ResonanceSlepton : public SUSYResonanceWidths {

public:

  // Constructor.
  ResonanceSlepton(int idResIn) {
      initBasic(idResIn);
}

private:

  bool getChannels(int idPDG) override;

  // Locally stored properties and couplings.

  // Initialize constants.
  virtual void initConstants();

  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Calculate width for currently considered channel.
  virtual void calcWidth(bool calledFromInit = false);

  double gf2;
  int id;

  // Three-body stau decaywidth classes
  One stauWidths;

};

}
