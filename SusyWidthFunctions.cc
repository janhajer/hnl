#include "SusyWidthFunctions.h"
#include "generic.hh"

namespace neutrino
{

// Function definitions (not found in the header) for the SUSY Resonance three-body decay width classes.
void ResonanceWidth::setPointers(Pythia8::ParticleData* particle_data_, Pythia8::CoupSM* standard_model_, Pythia8::Info* info_)
{
    particle_data = particle_data_;
    standard_model = standard_model_;
    info = info_;
}

// The StauWidths class. Width functions for 3-body stau decays.

double ResonanceWidth::getWidth(int idResIn, int to_id, int lepton_id)
{
    setChannel(idResIn, to_id, lepton_id);
// Calculate integration limits and return integrated width.
    if (idResIn % 2 == 0) return 0.;
    double width;
    return integrateGauss(width, 0., 1., 1e-3) ? width : 0.;
}

void ResonanceWidth::setChannel(int from_id, int to_id, int lepton_id)
{
    from_id = std::abs(from_id);
    to_id = std::abs(to_id);
    lepton_id = std::abs(lepton_id);

// Common masses.
    m_from = particle_data->m0(from_id);
    m_to = particle_data->m0(to_id);
    m_lepton = particle_data->m0(lepton_id);
    m_neutrino = particle_data->m0(9900014);

// Couplings etc.
    meson_lifetime = particle_data->mWidth(from_id);
    gf = standard_model->GF();

// Set function switch and internal propagators depending on decay product.
    if (to_id == 111) meson_channel = MesonChannel::Scalar;
    else if (to_id == 900111 || to_id == 113) meson_channel = MesonChannel::Vector;
    else if (to_id == 14 || to_id == 12) {
        m_to = particle_data->m0(to_id - 1);
        meson_channel = MesonChannel::None;
    } else {
        std::stringstream mess;
        mess << " unknown decay channel idIn = " << to_id;
        info->errorMsg("Warning in StauWidths::setChannel:", mess.str());
    }
    return;
}

double ResonanceWidth::f(double integrand)
{
// Decay width functions documented in arXiv:1212.2886 Citron et. al.
    double prefactor = meson_lifetime * neutrino_coupling * V2CKM * sqr(gf) / 64 / cube(M_PI) / sqr(m_from);
    double ffm = integrand;
    double ffp = integrand;

    double dalitz_min = (sqr(m_from) - sqr(m_lepton) - 2 * m_lepton * m_to - sqr(m_to) + sqr(m_neutrino)) / 2 / m_from;
    double dalitz_max = m_neutrino;
    double integral = dalitz_max - dalitz_min;
    double energy = (dalitz_min + dalitz_max) / 2;

    switch (meson_channel) {
    case MesonChannel::Scalar : {
        double term1 = sqr(ffm) * (integrand * (sqr(m_neutrino) + sqr(m_lepton)) - sqr(sqr(m_neutrino) - sqr(m_lepton)));
        double term2 = 2 * ffp * ffm * (sqr(m_neutrino) * (2 * sqr(m_from) - 2 * sqr(m_to) - 4 * energy * m_from - sqr(m_lepton) + sqr(m_neutrino) + integrand) + sqr(m_lepton) * (4 * energy * m_from + sqr(m_lepton) - sqr(m_neutrino) - integrand));
        double term3 = sqr(ffp) * ((4 * energy * m_from + sqr(m_lepton) - sqr(m_neutrino) - integrand) * (2 * sqr(m_from) - 2 * sqr(m_to) - 4 * energy * m_from - sqr(m_lepton) + sqr(m_neutrino) + integrand) - (2 * sqr(m_from) + 2 * sqr(m_to) - integrand) * (integrand - sqr(m_neutrino) - sqr(m_lepton)));
        return prefactor * integral * (term1 + term2 + term3);
    }
    case MesonChannel::Vector : {
        double term1 = sqrt((Pythia8::pow2(delta_mass) - qf2) * (Pythia8::pow2(delta_mass + 2 * m_neutrino) - qf2));
        double term2 = Pythia8::pow2(qf2 - Pythia8::pow2(m_to)) * (qf2 + Pythia8::pow2(m_to)) / (qf2 * qf2 * (Pythia8::pow2(qf2 - Pythia8::pow2(m_lepton)) + Pythia8::pow2(m_lepton * meson_lifetime)));
        return prefactor * (term1 * term2 * (term3 + term4));
    }
    default : {
        std::stringstream mess;
        mess << " unknown decay channel fnSwitch = " << fnSwitch;
        info->errorMsg("Warning in StauWidths::function:", mess.str());
        return 0.;
    }
    }
}



double One::getWidth(int idResIn, int to_id, int lepton_id)
{
    setChannel(idResIn, to_id, lepton_id);
// Calculate integration limits and return integrated width.
    if (idResIn % 2 == 0) return 0.;
    double width;
    return prefactor * (integrateGauss(width, sqr(y_l + y_N), sqr(1 - y_h), 1e-3 * y_h) ? width : 0.);
}


void One::setChannel(int from_id, int to_id, int lepton_id)
{
    from_id = std::abs(from_id);
    to_id = std::abs(to_id);
    lepton_id = std::abs(lepton_id);

// Common masses.
    m_from = particle_data->m0(from_id);
    m_to = particle_data->m0(to_id);
    m_lepton = particle_data->m0(lepton_id);
    m_neutrino = particle_data->m0(9900014);

    m_h = m_to;
    y_h = m_h / m_from;
    y_l = m_lepton / m_from;
    y_N = m_neutrino / m_from;

// Couplings etc.
    gf = standard_model->GF();

    prefactor = sqr(gf) * Pythia8::pow5(m_h) / 64 / cube(M_PI);

// Set function switch and internal propagators depending on decay product.
    channel = MesonChannel::Scalar;
}

double One::lambda(double a, double b, double c)
{
    return sqr(a) + sqr(b) + sqr(c) - 2 * a * b - 2 * a * c - 2 * b * c;
}

double One::Lambda(double xi)
{
    return std::sqrt(lambda(1, sqr(y_h), xi)) * std::sqrt(lambda(xi, sqr(y_N), sqr(y_l)));
}

double One::form_factor_plus(double sqr_q)
{
    return 1.;
}

double One::form_factor_zero(double sqr_q)
{
    return 1.;
}

double One::G_minus(double xi)
{
    return xi * (sqr(y_N) + sqr(y_l)) - sqr(sqr(y_N) - sqr(y_l));
}

double One::f(double integrand)
{
    switch (channel) {
    case MesonChannel::Scalar :{
        double term1 = sqr(form_factor_plus(integrand * sqr(m_h))) * cube(Lambda(integrand)) / 3 / cube(integrand);
        double term2 = sqr(form_factor_plus(integrand * sqr(m_h))) * Lambda(integrand) * G_minus(integrand) * lambda(1, sqr(y_h), integrand) / 2 / cube(integrand);
        double term3 = sqr(form_factor_zero(integrand * sqr(m_h))) * Lambda(integrand) * G_minus(integrand) * sqr(1 - sqr(y_h)) / 2 / cube(integrand);
        return term1 + term2 + term3;
    }
    case MesonChannel::Vector : return sqr(m_h) * sqr(y_h) * sqr(form_factor_zero(integrand * sqr(m_h))) * Lambda(integrand) * G_minus(integrand) * sqr(1 - sqr(y_h)) / 2 / cube(integrand);
    default : std::stringstream mess;
//         mess << " unknown decay channel fnSwitch = " << counter;
        info->errorMsg("Warning in StauWidths::function:", mess.str());
        return 0.;
    }
}

}
