#include "ThreeBodyWidth.hh"
#include "ThreeBodyWidth.hh"

namespace neutrino
{

namespace
{

const bool debug = false;

template<typename Object>
auto sqr(Object const& object) noexcept
{
    return object * object;
}

template<typename Object>
auto cube(Object const& object) noexcept
{
    return object * object * object;
}

void print() noexcept
{
    std::cout << std::endl;
}

template<typename Object, typename ... Arguments>
void print(Object const& object, Arguments ... arguments) noexcept
{
    std::cout << std::boolalpha << object << ' ';
    print(arguments ...);
}

bool is_vector(int id)
{
    return id == 113 || id == 213 || id == 313 || id == 323 || id == 413 || id == 423 || id == 433;
}

bool is_B_c_star(int id)
{
    return id == 543;
}

bool is_B_c(int id)
{
    return id == 541;
}

bool is_B_s_star(int id)
{
    return id == 533;
}

bool is_B_s(int id)
{
    return id == 531;
}

bool is_B_star(int id)
{
    return id == 513 || id == 523;
}

bool is_B(int id)
{
    return id == 511 || id == 521;
}

bool is_D_s_star(int id)
{
    return id == 433;
}

bool is_D_s(int id)
{
    return id == 431;
}

bool is_D_star(int id)
{
    return id == 413 || id == 423;
}

bool is_D(int id)
{
    return id == 411 || id == 421;
}

bool is_K_star(int id)
{
    return id == 313 || id == 323 || id == 311 || id == 321;
}

bool is_K(int id)
{
    return id == 130 || id == 310 || id == 311 || id == 321;
}

bool is_rho(int id)
{
    return id == 113 || id == 213;
}

bool is_pi(int id)
{
    return id == 111 || id == 211;
}

bool is_eta(int id)
{
    return id == 221 || id == 331;
}

bool is_neutral_Kaon(int id)
{
    return id == 130 || id == 310 || id == 311;
}

}

// void ThreeBody::set_pointers(Pythia8::ParticleData* particle_data_)
// {
//     if (debug) print("set pointers");
//     particle_data = particle_data_;
// }
//
// double ThreeBody::f(std::vector<double> integrands)
// {
//     return function(integrands[0]);
// }
//
// bool ThreeBody::integrate(double& result, double from, double to, double tolerance)
// {
//     std::vector<double> args(1);
//     return integrateGauss(result, 0, from, to, args, tolerance);
// }

void ThreeBodyWidth::set_pointers(Pythia8::ParticleData* particle_data_)
{
    if (debug) print("set pointers");
    particle_data = particle_data_;
}

double ThreeBodyWidth::get_width(int from_id_, int neutrino_id, int to_id_, int lepton_id_)
{
    if (debug) print("get_width", from_id_, "to", to_id_, "with", lepton_id_, "and", neutrino_id);

    id_from = std::abs(from_id_);
    id_to = std::abs(to_id_);
    auto id_lepton = std::abs(lepton_id_);
    auto id_neutrino = std::abs(neutrino_id);

    if (id_from < 100 || id_from > 1000) print("id from", id_from);
    if (id_to < 100 || id_to > 1000) print("id to", id_to);
    if (id_lepton < 10 || id_lepton > 20) print("id lepton", id_lepton);
    if (id_neutrino < 1E5) print("id neutrino", id_neutrino);

    if (id_to > id_from && id_from != 130 && id_to != 211) print("decaying", id_from, "to", id_to);
    if (id_lepton > 20 || id_lepton < 10) print("lepton?", id_lepton);

    mHat = particle_data->m0(id_from);
    m_to = particle_data->m0(id_to);

    auto y_h = m_to / mHat;
    auto y_l = particle_data->m0(id_lepton) / mHat;
    auto y_N = particle_data->m0(id_neutrino) / mHat;

    mr_h = sqr(y_h);
    mr_l = sqr(y_l);
    mr_N = sqr(y_N);

    double width;
    return integrate(width, sqr(y_l + y_N), sqr(1. - y_h), 1e-3 * y_N) ? width : 0.;
}

double ThreeBodyWidth::f(std::vector<double> integrands)
{
    return function(integrands[0]);
}

bool ThreeBodyWidth::integrate(double& result, double from, double to, double tolerance)
{
    std::vector<double> args(1);
    return integrateGauss(result, 0, from, to, args, tolerance);
}

double lambda(int id_from, int id_to, bool charged)
{
    if (is_neutral_Kaon(id_from) && id_to == 211) return charged ? 0.0267 : 0.0117;
    if (id_from == 321 && id_to == 111) return charged ? 0.0277 : 0.0183;
    print("not a valid kaon decay", id_from, id_to, charged);
    return 0;
}

double ThreeBodyWidth::K_form_factor(double q2, bool charged)
{
    return 0.970 * (1 + neutrino::lambda(id_from, id_to, charged) * q2 / sqr(particle_data->m0(211)));
}

double FF_d_0(int id)
{
    if (is_pi(id)) return 0.6117;
    if (is_K(id)) return 0.7647;
    print("not a proper D decay", id);
    return 0;
}

double ThreeBodyWidth::z(double q2)
{
    double t_p = sqr(mHat + m_to);
    double t_0 = (mHat + m_to) * sqr(std::sqrt(mHat) - std::sqrt(m_to));
    double first = std::sqrt(t_p - q2);
    double second = std::sqrt(t_p - t_0);
    return (first - second) / (first + second);
}

double c(int id, bool charged)
{
    if (is_pi(id)) return charged ? 1.985 : 1.188;
    if (is_K(id)) return charged ? 0.066 : 2.084;
    print("not a proper D decay", id, charged);
    return 0;
}

double P(int id, bool charged) // GeV^-2
{
    if (is_pi(id)) return charged ? 0.1314 : 0.0342;
    if (is_K(id)) return charged ? 0.224 : 0;
    print("not a proper D decay", id, charged);
    return 0;
}

double D_s_form_factor(double q2, bool charged)
{
    double f_p_eta = 0.495;
    double alpha_p_eta = 0.198;
    double m_D_s = 2.112; // GeV
    double f_0_eta = f_p_eta;
    double alpha_0_eta = 0;
    double q_r = q2 / sqr(m_D_s);
    return charged ? f_p_eta / (1 - q_r) / (1 - alpha_p_eta * q_r) :  f_0_eta / (1 - alpha_0_eta * q_r);
}

double ThreeBodyWidth::D_form_factor(double q2, bool charged)
{
    return (FF_d_0(id_to) - c(id_to, charged) * (z(q2) - z(0)) * (1. + (z(q2) + z(0)) / 2.)) / (1. - P(id_to, charged) * q2);
}

double m_pole_prefactor(double q2, int id, bool charged)
{
    if (is_D(id) || is_D_s(id)) return 1.;
    if (is_pi(id) || is_K(id)) {
        double xi = q2 / sqr(charged ? 5.325 : 5.65);
        return 1. / (1. - xi);
    }
    print("not a proper B decay", id, charged);
    return 0.;
}

double a_0(int id, bool charged)
{
    if (is_pi(id)) return charged ? 0.404 : 0.490;
    if (is_K(id)) return charged ? 0.360 : 0.233;
    if (is_D(id) || is_D_s(id)) return charged ? 0.909 : 0.794;
    print("not a proper B decay", id, charged);
    return 0.;
}

double a_1(int id, bool charged)
{
    if (is_pi(id)) return charged ? -0.68 : -1.61;
    if (is_K(id)) return charged ? -0.828 : 0.197;
    if (is_D(id) || is_D_s(id)) return charged ? -7.11 : -2.45;
    print("not a proper B decay", id, charged);
    return 0.;
}

double a_2(int id, bool charged)
{
    if (is_pi(id)) return charged ? -0.86 : 0.93;
    if (is_K(id)) return charged ? 1.1 : 0.18;
    if (is_D(id) || is_D_s(id)) return charged ? 66 : 33;
    print("not a proper B decay", id, charged);
    return 0.;
}

double a(int n, int id, int charged)
{
    switch (n) {
    case 0 : return a_0(id, charged);
    case 1 : return a_1(id, charged);
    case 2 : return a_2(id, charged);
    default : print("Case not covered");
    }
    return 0;
}

double zq2n(int n, double zq2, double zq2N3)
{
    switch (n) {
    case 0 : return 1.;
    case 1 : return zq2 - zq2N3;
    case 2 : return sqr(zq2) + 2 * zq2N3;
    default : print("unexpected integer");
    }
    return 0.;
}

double ThreeBodyWidth::B_form_factor(double q2, bool charged)
{
    double zq2 = z(q2);
    int N = 2;
    double zq2N3 = cube(zq2) / N;
    double sum = 0.;
    for (int n = 0; n < N; ++n) sum += a(n, id_to, charged) * zq2n(n, zq2, zq2N3);
    return m_pole_prefactor(q2, id_to, charged) * sum;
}

double ThreeBodyWidth::form_factor(double q2, bool charged)
{
    if (is_K(id_from)) return K_form_factor(q2, charged);
    if (is_D(id_from)) return D_form_factor(q2, charged);
    if (is_D_s(id_from)) return D_s_form_factor(q2, charged);
    if (is_B(id_from) || is_B_s(id_from)) return B_form_factor(q2, charged);
    print("not a proper pseudoscalar decay", id_from, id_to, charged, q2);
    return 0.;
}

double ThreeBodyWidth::form_factor_plus(double sqr_q)
{
    return form_factor(sqr_q, true);
}

double ThreeBodyWidth::form_factor_0(double sqr_q)
{
    return form_factor(sqr_q, false);
}

double fV(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 1.03;
    if (is_B(id_from) && is_D_star(id_to)) return 0.76;
    if (is_B(id_from) && is_rho(id_to)) return 0.295;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.95;
    if (is_B_s(id_from) && is_K_star(id_to)) return 0.291;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double fA0(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 0.76;
    if (is_B(id_from) && is_D_star(id_to)) return 0.69;
    if (is_B(id_from) && is_rho(id_to)) return 0.231;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.67;
    if (is_B_s(id_from) && is_K_star(id_to)) return 0.289;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double fA1(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 0.66;
    if (is_B(id_from) && is_D_star(id_to)) return 0.66;
    if (is_B(id_from) && is_rho(id_to)) return 0.269;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.70;
    if (is_B_s(id_from) && is_K_star(id_to)) return 0.287;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double fA2(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 0.49;
    if (is_B(id_from) && is_D_star(id_to)) return 0.62;
    if (is_B(id_from) && is_rho(id_to)) return 0.282;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.75;
    if (is_B_s(id_from) && is_K_star(id_to)) return 0.286;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double sigmaV(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 0.27;
    if (is_B(id_from) && is_D_star(id_to)) return 0.57;
    if (is_B(id_from) && is_rho(id_to)) return 0.875;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.372;
    if (is_B_s(id_from) && is_K_star(id_to)) return -0.516;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double sigmaA0(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 0.17;
    if (is_B(id_from) && is_D_star(id_to)) return 0.59;
    if (is_B(id_from) && is_rho(id_to)) return 0.796;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.350;
    if (is_B_s(id_from) && is_K_star(id_to)) return -0.383;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double sigmaA1(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 0.30;
    if (is_B(id_from) && is_D_star(id_to)) return 0.78;
    if (is_B(id_from) && is_rho(id_to)) return 0.54;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.463;
    if (is_B_s(id_from) && is_K_star(id_to)) return 0.;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double sigmaA2(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 0.67;
    if (is_B(id_from) && is_D_star(id_to)) return 1.40;
    if (is_B(id_from) && is_rho(id_to)) return 1.34;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 1.04;
    if (is_B_s(id_from) && is_K_star(id_to)) return 1.05;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double xiV(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 0.;
    if (is_B(id_from) && is_D_star(id_to)) return 0.;
    if (is_B(id_from) && is_rho(id_to)) return 0.;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.561;
    if (is_B_s(id_from) && is_K_star(id_to)) return 2.10;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double xiA0(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 0.;
    if (is_B(id_from) && is_D_star(id_to)) return 0.;
    if (is_B(id_from) && is_rho(id_to)) return 0.055;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.600;
    if (is_B_s(id_from) && is_K_star(id_to)) return 1.58;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double xiA1(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 0.20;
    if (is_B(id_from) && is_D_star(id_to)) return 0.;
    if (is_B(id_from) && is_rho(id_to)) return 0.;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.510;
    if (is_B_s(id_from) && is_K_star(id_to)) return 1.06;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double xiA2(int id_from, int id_to)
{
    if (is_D(id_from) && is_K_star(id_to)) return 0.16;
    if (is_B(id_from) && is_D_star(id_to)) return 0.41;
    if (is_B(id_from) && is_rho(id_to)) return -0.21;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 0.070;
    if (is_B_s(id_from) && is_K_star(id_to)) return -0.074;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double MP(int id_from, int id_to) // GeV
{
    if (is_D(id_from) && is_K_star(id_to)) return 1.969;
    if (is_B(id_from) && is_D_star(id_to)) return 6.275;
    if (is_B(id_from) && is_rho(id_to)) return 5.279;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 6.275;
    if (is_B_s(id_from) && is_K_star(id_to)) return 5.367;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double MV(int id_from, int id_to) // GeV
{
    if (is_D(id_from) && is_K_star(id_to)) return 2.112;
    if (is_B(id_from) && is_D_star(id_to)) return 6.331;
    if (is_B(id_from) && is_rho(id_to)) return 5.325;
    if (is_B_s(id_from) && is_D_s_star(id_to)) return 6.331;
    if (is_B_s(id_from) && is_K_star(id_to)) return 5.415;
    print("meson", id_from, id_to, "not found");
    return 0;
}

double V(double q2, int id_from, int id_to)
{
    double ratio = q2 / sqr(MV(id_from, id_to));
    return fV(id_from, id_to) / (1. - ratio) / (1. - sigmaV(id_from, id_to) * ratio - xiV(id_from, id_to) * sqr(ratio));
}

double A0(double q2, int id_from, int id_to)
{
    double ratioP = q2 / sqr(MP(id_from, id_to));
    double ratioV = q2 / sqr(MV(id_from, id_to));
    return fA0(id_from, id_to) / (1. - ratioP) / (1. - sigmaA0(id_from, id_to) * ratioV - xiA0(id_from, id_to) * sqr(ratioV));
}

double A1(double q2, int id_from, int id_to)
{
    double ratio = q2 / sqr(MV(id_from, id_to));
    return fA1(id_from, id_to) / (1. - sigmaA1(id_from, id_to) * ratio - xiA1(id_from, id_to) * sqr(ratio));
}

double A2(double q2, int id_from, int id_to)
{
    double ratio = q2 / sqr(MV(id_from, id_to));
    return fA2(id_from, id_to) / (1. - sigmaA2(id_from, id_to) * ratio - xiA2(id_from, id_to) * sqr(ratio));
}

double ThreeBodyWidth::g(double q2)
{
    return V(q2, id_from, id_to) / (mHat + m_to);
}

double ThreeBodyWidth::am(double q2)
{
    return (A2(q2, id_from, id_to) * (mHat - m_to) - A1(q2, id_from, id_to) * (mHat - m_to) + 2 * A0(q2, id_from, id_to) * m_to) / q2;
}

double ThreeBodyWidth::ff(double q2)
{
    return A1(q2, id_from, id_to) * (mHat + m_to);
}

double ThreeBodyWidth::ap(double q2)
{
    return - A2(q2, id_from, id_to) / (mHat + m_to);
}

namespace{

double lambda(double a, double b, double c)
{
    return sqr(a) + sqr(b) + sqr(c) - 2 * a * b - 2 * a * c - 2 * b * c;
}

}

double Lambda(double x, double mr_to, double mr_1, double mr_2)
{
    return std::sqrt(lambda(1., mr_to, x) * lambda(x, mr_1, mr_2));
}

double ThreeBodyWidth::Lambda(double xi)
{
    return neutrino::Lambda(xi, mr_h, mr_N, mr_l);
}

double ThreeBodyWidth::Gm(double xi)
{
    return xi * (mr_N + mr_l) - sqr(mr_N - mr_l);
}

double ThreeBodyWidth::Gp(double xi)
{
    return xi * (mr_N + mr_l) + sqr(mr_N - mr_l);
}

double ThreeBodyWidth::F(double xi)
{
    return sqr(1 - xi) - 2 * mr_h * (1 + xi) + sqr(mr_h);
}

double ThreeBodyWidth::function(double xi)
{
    double mHat2 = sqr(mHat);
    double q2 = xi * mHat2;
    double Lambdaxi = Lambda(xi);
    double Gmxi = Gm(xi);
    double xi2 = sqr(xi);
    double xi3 = xi2 * xi;
    if (!is_vector(id_to)) {
        double ffp2q2 = sqr(form_factor_plus(q2));
        double term1 = ffp2q2 * cube(Lambdaxi) / 3. / xi3;
        double term2 = ffp2q2 * Lambdaxi * Gmxi * lambda(1., mr_h, xi) / 2. / xi3;
        double term3 = sqr(form_factor_0(q2)) * Lambdaxi * Gmxi * sqr(1. - mr_h) / 2. / xi3;
        return term1 + term2 + term3;
    } else {
        double fq2 = ff(q2);
        double Fxi = F(xi);
        double Gpxi = Gp(xi);
        double apq2 = ap(q2);
        double amq2 = am(q2);
        double term1 = mHat2 * mr_h / 3. / xi2 * g(q2) * Lambdaxi * Fxi * (2 * xi2 - Gpxi);
        double term2 = 1. / 24. / mHat2 / xi3 * sqr(fq2) * Lambdaxi * (3. * Fxi * (xi2 - sqr(mr_l - mr_N)) - sqr(Lambdaxi) + 12. * mr_h * xi * (2 * xi2 - Gpxi));
        double term3 = mHat2 / 24. / xi3 * sqr(apq2) * Lambdaxi * Fxi * (Fxi * (2 * xi2 - Gpxi) + 3. * Gmxi * sqr(1. - mr_h));
        double term4 = mHat2 / 8. / xi * sqr(amq2) * Lambdaxi * Fxi * Gmxi;
        double term5 = 1. / 12. / xi3 * fq2 * apq2 * Lambdaxi * (3 * xi * Fxi * Gmxi + (1 - xi - mr_h) * (3 * Fxi * (xi2 - sqr(mr_l - mr_N)) - sqr(Lambdaxi)));
        double term6 = 1. / 4. / xi2 / fq2 * amq2 * Lambdaxi * Fxi * Gmxi;
        double term7 = mHat2 / 4. / xi2 * apq2 * amq2 * Lambdaxi * Fxi * Gmxi * (1. - mr_h);
        return term1 + term2 + term3 + term4 + term5 + term6 + term7;
    }
}


////////////////////////////







void NeutrinoThreeBodyWidth::set_pointers(Pythia8::ParticleData* particle_data_)
{
    if (debug) print("set pointers");
    particle_data = particle_data_;
}

double NeutrinoThreeBodyWidth::f(std::vector<double> integrands)
{
    return function(integrands[0]);
}

bool NeutrinoThreeBodyWidth::integrate(double& result, double from, double to, double tolerance)
{
    std::vector<double> args(1);
    return integrateGauss(result, 0, from, to, args, tolerance);
}

double NeutrinoThreeBodyWidth::get_width(int id_from, int id_lepton, int id_up, int id_down)
{
    if (debug) print("get_width", id_from, "to", id_up, "with", id_down, "and", id_lepton);

    if (id_from < 1E5) print("id neutrino", id_lepton);
    if (id_lepton < 10 || id_lepton > 20) print("id lepton", id_lepton);
    if (id_down > 20) print("id down", id_down);
    if (id_up > 20) print("id up", id_up);

    auto m_from = particle_data->m0(id_from);

    auto y_l = particle_data->m0(id_lepton) / m_from;
    auto y_u = particle_data->m0(id_up) / m_from;
    auto y_d = particle_data->m0(id_down) / m_from;

    mr1 = sqr(y_l);
    mr2 = sqr(y_d);
    mr3 = sqr(y_u);

    double width;
    return integrate(width, sqr(y_d + y_u), sqr(1. - y_l), 1e-3 * y_l) ? width : 0.;
}

double NeutrinoThreeBodyWidth::function(double x)
{
    return 12. / x * (x - mr2 - mr3) * (1 + mr1 - x) * Lambda(x, mr1, mr2, mr3);
}

}
