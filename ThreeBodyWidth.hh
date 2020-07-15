#pragma once

#include "Pythia8/ParticleData.h"
#include "Pythia8/StandardModel.h"

namespace neutrino
{

struct ThreeBodyWidth : public Pythia8::FunctionEncapsulator
{
    void set_pointers(Pythia8::ParticleData* particleDataPtrIn);
    double get_width(int from_id_, int to_id_, int lepton_id_);
protected:
    virtual double f(std::vector<double> integrands) override;
private:
    bool integrate(double& result, double from, double to, double tolerance);
    double function(double integrand);
    double lambda(double a, double b, double c);
    double Lambda(double xi);
    double ffp(double q2);
    double ff0(double q2);
    double Gm(double xi);
    double form_factor(double q2, int id, int id_2, bool charged);
    double K_form_factor(double q2, bool charged);
    double D_form_factor(double q2, bool charged);
    double B_form_factor(double q2, int id, bool charged);
    double F(double xi);
    double Gp(double xi);
    double z(double q2);
    double g(double q2);
    double am(double q2);
    double ff(double q2);
    double ap(double q2);
private:
    Pythia8::ParticleData* particle_data;
    double m_from;
    double m_to;
    double mr_h;
    double mr_N;
    double mr_l;
    int id_from;
    int id_to;
};

}
