#pragma once

#include "Pythia8/ParticleData.h"
#include "Pythia8/StandardModel.h"

namespace neutrino
{

enum class MesonChannel {None, Scalar, Vector};

struct ThreeBodyWidth : public Pythia8::FunctionEncapsulator
{
    void setPointers(Pythia8::ParticleData* particleDataPtrIn);
    double getWidth(int from_id_, int to_id_, int lepton_id_);
protected:
    virtual double f(std::vector<double> integrands) override;
private:
    MesonChannel meson_channel;
    bool integrate(double& result, double from, double to, double tolerance);
    Pythia8::ParticleData* particle_data;
    int id_from;
    int id_to;
    int id_lepton;
    double m_from;
    double m_lepton;
    double m_neutrino;
    double m_to;

    void setChannel(int from_id, int to_id, int lepton_id) ;
    double function(double integrand);
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
