// SusyResonanceWidths.h is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand
// Main author of this file: N. Desai
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for resonance properties: dynamical widths etc.
// WidthFunctions: base class for SUSY 3-body decay width functions.

#pragma once

#include "Pythia8/ParticleData.h"
#include "Pythia8/StandardModel.h"
#include <complex>

namespace neutrino {

enum class MesonChannel {None, Scalar, Vector};

class ResonanceWidth : public Pythia8::FunctionEncapsulator {
public:
  ResonanceWidth(double neutrino_coupling_) :
  particle_data(),
  standard_model(),
  info(),
  from_id(),
  idInt(),
  id1(),
  id2(),
  id3(),
  id4(),
  m_from(),
  m_lepton(),
  m_neutrino(),
  m_to(),
  meson_lifetime(),
  neutrino_coupling(neutrino_coupling_)
  {};
  void setPointers( Pythia8::ParticleData* particleDataPtrIn, Pythia8::CoupSM* coupSUSYPtrIn, Pythia8::Info* infoPtrIn);
  double getWidth(int idResIn, int idIn);
protected:
  // Wrappers to simplify using FunctionEncapsulator's integrator.
  void setChannel(int idResIn, int idIn, int lepton_id) ;
  double f(double integrand) ;
  double f(std::vector<double> integrands) override { return f(integrands[0]); }
  bool integrateGauss(double& result, double xLo, double xHi, double tol)  {
    std::vector<double> tmp(1);
    return FunctionEncapsulator::integrateGauss(result, 0, xLo, xHi, tmp, tol);
  }
  Pythia8::ParticleData* particle_data;
  Pythia8::CoupSM* standard_model;
  Pythia8::Info* info;
  int from_id;
  int idInt;
  int id1;
  int id2;
  int id3;
  int id4;
  double m_from;
  double m_lepton;
  double m_neutrino;
  double m_to;
  double meson_lifetime;
  double neutrino_coupling;
  MesonChannel meson_channel; // Switch between multiple functions
  double delta_mass;
  double gf;
  double cons;
  double wparam;
  std::complex<double> gL;
  std::complex<double> gR;
};


struct One : ResonanceWidth {
  void setChannel(int from_id, int to_id, int lepton_id) ;
  double getWidth(int from_id, int to_id);
  double f(double integrand);
  double lambda(double a, double b, double c);
  double Lambda(double xi);
  double form_factor_plus(double sqr_q);
  double form_factor_zero(double sqr_q);
  double G_minus(double xi);
  MesonChannel channel;
  double m_h;
  double y_h;
  double y_N;
  double y_l;
  double prefactor;
};

}
