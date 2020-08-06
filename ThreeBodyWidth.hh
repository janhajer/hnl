#pragma once

#include <vector>

namespace Pythia8{
    class ParticleData;
    class Settings;
}

namespace hnl
{

struct ThreeBodyWidth
{
    virtual ~ThreeBodyWidth() = default;
protected:
    virtual double function(std::vector<double> const& integrands) const;
    virtual double function(double integrand) const = 0;
    virtual bool integrate(double& result, double from, double to, double tolerance) const;
    bool integrateGauss(double& result, int iArg, double xLo, double xHi, std::vector<double> args, double tol = 1e-6) const;
    bool brent(double& solution, double target, int iArg, double xLo, double xHi, std::vector<double> argsIn, double tol = 1e-6, int maxIter = 10000) const;
};

struct MesonThreeBodyWidth : public ThreeBodyWidth {
    void set_pointers(Pythia8::ParticleData* particleDataPtrIn);
    double get_width(int from_id_, int neutrino_id, int to_id_, int lepton_id_);
protected:
    double function(double integrand) const override;
private:
    double Lambda(double xi) const;
    double F(double xi) const;
    double Gp(double xi) const;
    double Gm(double xi) const;
    double form_factor_plus(double q2) const;
    double form_factor_0(double q2) const;
    double form_factor(double q2, bool charged) const;
    double K_form_factor(double q2, bool charged) const;
    double D_form_factor(double q2, bool charged) const;
    double B_form_factor(double q2, bool charged) const;
    double z(double q2) const;
    double g(double q2) const;
    double am(double q2) const;
    double ff(double q2) const;
    double ap(double q2) const;
private:
    Pythia8::ParticleData* particle_data;
    int id_from;
    int id_to;
    double mHat;
    double m_to;
    double mr_h;
    double mr_l;
    double mr_N;
};

struct NeutrinoThreeBodyWidth : public ThreeBodyWidth {
    void set_pointers(Pythia8::ParticleData* particle_data_, Pythia8::Settings* settings_);
    double get_width(int from_id_, int neutrino_id, int to_id_, int lepton_id_);
protected:
    double function(double integrand) const override;
private:
    double get_mass(int id) const;
private:
    Pythia8::ParticleData* particle_data;
    Pythia8::Settings* settings;
    double mr1;
    double mr2;
    double mr3;
};

}
