#pragma once

#include <vector>

namespace Pythia8{
    class ParticleData;
}

namespace neutrino
{

// Base class to encapsulate a (double) function of an arbitrary number
// of (double) arguments (to avoid using function pointers).
class FunctionEncapsulator
{
public:
    FunctionEncapsulator() {};
    virtual ~FunctionEncapsulator() { };
    virtual double f(std::vector<double> args);
    // Integrate over function argument iArg, using Gaussian quadrature.
    bool integrateGauss(double& result, int iArg, double xLo, double xHi, std::vector<double> args, double tol = 1e-6);
    // Solve f(args) = target for argument iArg, using Brent's method
    bool brent(double& solution, double target, int iArg, double xLo, double xHi, std::vector<double> argsIn, double tol = 1e-6, int maxIter = 10000);
};

// struct ThreeBody : public Pythia8::FunctionEncapsulator
// {
//     void set_pointers(Pythia8::ParticleData* particleDataPtrIn);
// protected:
//     virtual double f(std::vector<double> integrands) override;
//     virtual bool integrate(double& result, double from, double to, double tolerance);
//     virtual double function(double integrand) = 0;
//     Pythia8::ParticleData* particle_data;
//     double mr1;
//     double mr2;
//     double mr3;
// };

struct ThreeBodyWidth : public FunctionEncapsulator {
    void set_pointers(Pythia8::ParticleData* particleDataPtrIn);
    double get_width(int from_id_, int neutrino_id, int to_id_, int lepton_id_);
protected:
    virtual double f(std::vector<double> integrands) override;
private:
    bool integrate(double& result, double from, double to, double tolerance);
    double function(double integrand);
    double Lambda(double xi);
    double F(double xi);
    double Gp(double xi);
    double Gm(double xi);
    double form_factor_plus(double q2);
    double form_factor_0(double q2);
    double form_factor(double q2, bool charged);
    double K_form_factor(double q2, bool charged);
    double D_form_factor(double q2, bool charged);
    double B_form_factor(double q2, bool charged);
    double z(double q2);
    double g(double q2);
    double am(double q2);
    double ff(double q2);
    double ap(double q2);
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

struct NeutrinoThreeBodyWidth : public FunctionEncapsulator {
    void set_pointers(Pythia8::ParticleData* particleDataPtrIn);
    double get_width(int from_id_, int neutrino_id, int to_id_, int lepton_id_);
protected:
    virtual double f(std::vector<double> integrands) override;
private:
    bool integrate(double& result, double from, double to, double tolerance);
    double function(double integrand);
private:
    Pythia8::ParticleData* particle_data;
    double mr1;
    double mr2;
    double mr3;
};

}
