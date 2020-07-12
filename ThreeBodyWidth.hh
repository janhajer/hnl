#pragma once

#include "Pythia8/ParticleData.h"
#include "Pythia8/StandardModel.h"

namespace neutrino
{

enum class MesonChannel {None, Scalar, Vector};

class ThreeBodyWidth : public Pythia8::FunctionEncapsulator
{
public:
    void setPointers(Pythia8::ParticleData* particleDataPtrIn, Pythia8::CoupSM* coupSUSYPtrIn, Pythia8::Info* infoPtrIn);
    double getWidth(int from_id_, int to_id_, int lepton_id_);
    MesonChannel meson_channel;
protected:
    double f(std::vector<double> integrands) override
    {
        return f(integrands[0]);
    }
    bool integrateGauss(double& result, double xLo, double xHi, double tol)
    {
        std::vector<double> tmp(1);
        return FunctionEncapsulator::integrateGauss(result, 0, xLo, xHi, tmp, tol);
    }
    Pythia8::ParticleData* particle_data;
    Pythia8::CoupSM* standard_model;
    Pythia8::Info* info;
    int id_from;
    int id_to;
    int id_lepton;
    double m_from;
    double m_lepton;
    double m_neutrino;
    double m_to;

    void setChannel(int from_id, int to_id, int lepton_id) ;
    double f(double integrand);
    double lambda(double a, double b, double c);
    double Lambda(double xi);
    double ffp(double sqr_q);
    double ff0(double sqr_q);
    double Gm(double xi);
    double form_factor(double q_sqr, int id, int id_2, bool charged);
    double K_meson_form_factor(double q_sqr, int id, bool charged);
    double D_meson_form_factor(double q_sqr, int id, bool charged);
    double D_meson_eta_form_factor(double q_sqr, bool charged);
    double B_meson_form_factor(double q_sqr, int id, bool charged);
    double F(double xi);
    double Gp(double xi);
    MesonChannel channel;
    double y_h;
    double y_N;
    double y_l;
    double mr_h;
    double mr_N;
    double mr_l;
    double prefactor;
    double z(double q_sqr);
    double g(double q2);
    double am(double q2);
    double ff(double q2);
    double ap(double q2);
};

}
