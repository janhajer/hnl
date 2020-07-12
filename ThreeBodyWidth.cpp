#include "ThreeBodyWidth.hh"
// #include "generic.hh"


namespace neutrino
{

namespace
{

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

}

// Function definitions (not found in the header) for the SUSY Resonance three-body decay width classes.
void ThreeBodyWidth::setPointers(Pythia8::ParticleData* particle_data_, Pythia8::CoupSM* standard_model_, Pythia8::Info* info_)
{
    print("set pointers");
    particle_data = particle_data_;
    standard_model = standard_model_;
    info = info_;
}

double ThreeBodyWidth::getWidth(int from_id_, int to_id_, int lepton_id_)
{
    print("get_width");
    setChannel(from_id_, to_id_, lepton_id_);
    double width;
    return integrateGauss(width, sqr(y_l + y_N), sqr(1 - y_h), 1e-3 * y_N) ? width : 0.;
}

void ThreeBodyWidth::setChannel(int from_id_, int to_id_, int lepton_id_)
{
    print("set channel from ", from_id_, "to", to_id_, "with", lepton_id_);
    id_from = std::abs(from_id_);
    id_to = std::abs(to_id_);
    id_lepton = std::abs(lepton_id_);

    m_from = particle_data->m0(id_from);
    m_to = particle_data->m0(id_to);
    m_lepton = particle_data->m0(id_lepton);
    m_neutrino = particle_data->m0(9900014);

    y_h = m_to / m_from;
    y_l = m_lepton / m_from;
    y_N = m_neutrino / m_from;

    mr_h = sqr(y_h);
    mr_l = sqr(y_l);
    mr_N = sqr(y_N);

    channel = id_to == 113 || id_to == 213 || id_to == 313 || id_to == 323 || id_to == 413 || id_to == 423 || id_to == 433 ? MesonChannel::Vector : MesonChannel::Scalar;
}

double ThreeBodyWidth::lambda(double a, double b, double c)
{
    return sqr(a) + sqr(b) + sqr(c) - 2 * a * b - 2 * a * c - 2 * b * c;
}

double ThreeBodyWidth::Lambda(double xi)
{
    return std::sqrt(lambda(1, mr_h, xi)) * std::sqrt(lambda(xi, mr_N, mr_l));
}

double lambda(int id, bool charged)
{
    return charged ?
           (id == 311 ? 0.0267 : 0.0277) :
           (id == 311 ? 0.0117 : 0.0183);
}

double ThreeBodyWidth::K_meson_form_factor(double q2, int id, bool charged)
{
    return 0.970 * (1 + neutrino::lambda(id, charged) * q2 / sqr(particle_data->m0(211)));
}

double FF_d_0(int id)
{
    return id < 300 ? 0.6114 : 0.7647;
}

double ThreeBodyWidth::z(double q2)
{
    double t_p = sqr(m_from + m_to);
    double t_0 = (m_from + m_to) * sqr(std::sqrt(m_from) - std::sqrt(m_to));
    double first = std::sqrt(t_p - q2);
    double second = std::sqrt(t_p - t_0);
    return (first - second) / (first + second);
}

double c(int id, bool charged)
{
    return charged ?
           (id > 300 ? 0.066 : 1.985) :
           (id > 300 ? 2.084 : 1.188);
}

double P(int id, bool charged)
{
    return charged ?
           (id > 300 ? 0.224 : 0.1314) :
           (id > 300 ? 0 : 0.0342);
}

double ThreeBodyWidth::D_meson_eta_form_factor(double q2, bool charged)
{
    double f_p_eta = 0.495;
    double alpha_p_eta = 0.198;
    double m_D_s = 2.112; // GeV
    double f_0_eta = f_p_eta;
    double alpha_0_eta = 0;
    double q_r = q2 / sqr(m_D_s);
    return charged ? f_p_eta / ((1 - q_r) * (1 - alpha_p_eta * q_r)) :  f_0_eta / (1 - alpha_0_eta * q_r);
}

double ThreeBodyWidth::D_meson_form_factor(double q2, int id, bool charged)
{
    if (id == 221 || id == 331) return D_meson_eta_form_factor(q2, charged);
    return (FF_d_0(id) * c(id, charged) * (z(q2) - z(0)) * (1 + z(q2) - z(0) / 2)) / (1 - P(id, charged) * q2);
}
double m_pole_prefactor(double q2, int id, bool charged)
{
    if (id < 400) return 1. / (1. - q2 / sqr(charged ? 5.325 : 5.65));
    return 1;
}

double a_0(int id, bool charged)
{
    if (id < 300) return charged ? 0.404 : 0.490;
    if (id < 400) return charged ? 0.360 : 0.233;
    return charged ? 0.909 : 0.794;
}

double a_1(int id, bool charged)
{
    if (id < 300) return charged ? -0.68 : -1.61;
    if (id < 400) return charged ? -0.828 : 0.197;
    return charged ? -7.11 : -2.45;
}

double a_2(int id, bool charged)
{
    if (id < 300) return charged ? -0.86 : 0.93;
    if (id < 400) return charged ? 1.1 : 0.18;
    return charged ? 66 : 33;
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

double ThreeBodyWidth::B_meson_form_factor(double q2, int id, bool charged)
{
    double z_q = z(q2);
    int N = 3;
    double sum = 0;
    for (auto n = 0; n < N; ++n) sum += a(n, id, charged) * (std::pow(z_q, n) - std::pow(-1, n - N) * n / N * std::pow(z_q, N));
    return m_pole_prefactor(q2, id, charged) * sum;
}

auto get_hundredth(int id)
{
    return id > 999 ? 100 * ((id / 100) % 10) + 10 * ((id / 10) % 10) + id % 10 : id;
}

double ThreeBodyWidth::form_factor(double q2, int from_id, int to_id, bool charged)
{
    int id = get_hundredth(from_id);
    if (id < 400) return K_meson_form_factor(q2, to_id, charged);
    if (id < 500) return D_meson_form_factor(q2, to_id, charged);
    return B_meson_form_factor(q2, to_id, charged);
}

double ThreeBodyWidth::ffp(double sqr_q)
{
    return form_factor(sqr_q, id_from, id_to, true);
}

double ThreeBodyWidth::ff0(double sqr_q)
{
    return form_factor(sqr_q, id_from, id_to, false);
}

double ThreeBodyWidth::Gm(double xi)
{
    return xi * (mr_N + mr_l) - sqr(mr_N - mr_l);
}

double ThreeBodyWidth::F(double xi)
{
    return sqr(1 - xi) - 2 * mr_h * (1 + xi) + sqr(mr_h);
}

double ThreeBodyWidth::Gp(double xi)
{
    return xi * (mr_N + mr_l) + sqr(mr_N - mr_l);
}

bool is_D(int id)
{
    id = get_hundredth(id);
    return id > 400 && id < 430;
}

bool is_B(int id)
{
    id = get_hundredth(id);
    return id > 500 && id < 530;
}

bool is_Bs(int id)
{
    id = get_hundredth(id);
    return id > 530 && id < 540;
}

bool is_Ds(int id)
{
    id = get_hundredth(id);
    return id > 430 && id < 440;
}

bool is_K(int id)
{
    id = get_hundredth(id);
    return id > 300 && id < 400;
}

bool is_rho(int id)
{
    id = get_hundredth(id);
    return id < 300;
}

double fV(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 1.03;
    if (is_B(id1) && is_D(id2)) return 0.76;
    if (is_B(id1) && is_rho(id2)) return 0.295;
    if (is_Bs(id1) && is_Ds(id2)) return 0.95;
    if (is_Bs(id1) && is_K(id2)) return 0.291;
    print("meson", id1, id2, "not found");
    return 0;
}

double fA0(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 0.76;
    if (is_B(id1) && is_D(id2)) return 0.69;
    if (is_B(id1) && is_rho(id2)) return 0.231;
    if (is_Bs(id1) && is_Ds(id2)) return 0.67;
    if (is_Bs(id1) && is_K(id2)) return 0.289;
    print("meson", id1, id2, "not found");
    return 0;
}

double fA1(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 0.66;
    if (is_B(id1) && is_D(id2)) return 0.66;
    if (is_B(id1) && is_rho(id2)) return 0.269;
    if (is_Bs(id1) && is_Ds(id2)) return 0.70;
    if (is_Bs(id1) && is_K(id2)) return 0.287;
    print("meson", id1, id2, "not found");
    return 0;
}

double fA2(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 0.49;
    if (is_B(id1) && is_D(id2)) return 0.62;
    if (is_B(id1) && is_rho(id2)) return 0.282;
    if (is_Bs(id1) && is_Ds(id2)) return 0.75;
    if (is_Bs(id1) && is_K(id2)) return 0.286;
    print("meson", id1, id2, "not found");
    return 0;
}

double sigmaV(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 0.27;
    if (is_B(id1) && is_D(id2)) return 0.57;
    if (is_B(id1) && is_rho(id2)) return 0.875;
    if (is_Bs(id1) && is_Ds(id2)) return 0.372;
    if (is_Bs(id1) && is_K(id2)) return -0.516;
    print("meson", id1, id2, "not found");
    return 0;
}

double sigmaA0(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 0.17;
    if (is_B(id1) && is_D(id2)) return 0.59;
    if (is_B(id1) && is_rho(id2)) return 0.796;
    if (is_Bs(id1) && is_Ds(id2)) return 0.350;
    if (is_Bs(id1) && is_K(id2)) return -0.383;
    print("meson", id1, id2, "not found");
    return 0;
}

double sigmaA1(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 0.30;
    if (is_B(id1) && is_D(id2)) return 0.78;
    if (is_B(id1) && is_rho(id2)) return 0.54;
    if (is_Bs(id1) && is_Ds(id2)) return 0.463;
    if (is_Bs(id1) && is_K(id2)) return 0.;
    print("meson", id1, id2, "not found");
    return 0;
}

double sigmaA2(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 0.67;
    if (is_B(id1) && is_D(id2)) return 1.40;
    if (is_B(id1) && is_rho(id2)) return 1.34;
    if (is_Bs(id1) && is_Ds(id2)) return 1.04;
    if (is_Bs(id1) && is_K(id2)) return 1.05;
    print("meson", id1, id2, "not found");
    return 0;
}

double xiV(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 0.;
    if (is_B(id1) && is_D(id2)) return 0.;
    if (is_B(id1) && is_rho(id2)) return 0.;
    if (is_Bs(id1) && is_Ds(id2)) return 0.561;
    if (is_Bs(id1) && is_K(id2)) return 2.10;
    print("meson", id1, id2, "not found");
    return 0;
}

double xiA0(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 0.;
    if (is_B(id1) && is_D(id2)) return 0.;
    if (is_B(id1) && is_rho(id2)) return 0.055;
    if (is_Bs(id1) && is_Ds(id2)) return 0.600;
    if (is_Bs(id1) && is_K(id2)) return 1.58;
    print("meson", id1, id2, "not found");
    return 0;
}

double xiA1(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 0.20;
    if (is_B(id1) && is_D(id2)) return 0.;
    if (is_B(id1) && is_rho(id2)) return 0.;
    if (is_Bs(id1) && is_Ds(id2)) return 0.510;
    if (is_Bs(id1) && is_K(id2)) return 1.06;
    print("meson", id1, id2, "not found");
    return 0;
}

double xiA2(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 0.16;
    if (is_B(id1) && is_D(id2)) return 0.41;
    if (is_B(id1) && is_rho(id2)) return -0.21;
    if (is_Bs(id1) && is_Ds(id2)) return 0.070;
    if (is_Bs(id1) && is_K(id2)) return -0.074;
    print("meson", id1, id2, "not found");
    return 0;
}

double MP(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 1.969;
    if (is_B(id1) && is_D(id2)) return 6.275;
    if (is_B(id1) && is_rho(id2)) return 5.279;
    if (is_Bs(id1) && is_Ds(id2)) return 6.275;
    if (is_Bs(id1) && is_K(id2)) return 5.367;
    print("meson", id1, id2, "not found");
    return 0;
}

double MV(int id1, int id2)
{
    if (is_D(id1) && is_K(id2)) return 2.112;
    if (is_B(id1) && is_D(id2)) return 6.331;
    if (is_B(id1) && is_rho(id2)) return 5.325;
    if (is_Bs(id1) && is_Ds(id2)) return 6.331;
    if (is_Bs(id1) && is_K(id2)) return 5.415;
    print("meson", id1, id2, "not found");
    return 0;
}

double V(double q2, int id1, int id2)
{
    double ratio = q2 / sqr(MV(id1, id2));
    return fV(id1, id2) / (1. - ratio) / (1. - sigmaV(id1, id2) * ratio - xiV(id1, id2) * sqr(ratio));
}

double A0(double q2, int id1, int id2)
{
    double ratioP = q2 / sqr(MP(id1, id2));
    double ratioV = q2 / sqr(MV(id1, id2));
    return fA0(id1, id2) / (1. - ratioP) / (1. - sigmaA0(id1, id2) * ratioV - xiA0(id1, id2) * sqr(ratioV));
}

double A1(double q2, int id1, int id2)
{
    double ratio = q2 / sqr(MV(id1, id2));
    return fA1(id1, id2) / (1. - sigmaA1(id1, id2) * ratio - xiA1(id1, id2) * sqr(ratio));
}

double A2(double q2, int id1, int id2)
{
    double ratio = q2 / sqr(MV(id1, id2));
    return fA2(id1, id2) / (1. - sigmaA2(id1, id2) * ratio - xiA2(id1, id2) * sqr(ratio));
}

double ThreeBodyWidth::g(double q2)
{
    return V(q2, id_from, id_to) / (m_from + m_to);
}

double ThreeBodyWidth::am(double q2)
{
    return (A2(q2, id_from, id_to) * (m_from - m_to) - A1(q2, id_from, id_to) * (m_from - m_to) + 2 * A0(q2, id_from, id_to) * m_to) / q2;
}

double ThreeBodyWidth::ff(double q2)
{
    return A1(q2, id_from, id_to) / (m_from + m_to);
}

double ThreeBodyWidth::ap(double q2)
{
    return - A2(q2, id_from, id_to) / (m_from + m_to);
}

double ThreeBodyWidth::f(double xi)
{
    double m_from_2 = sqr(m_from);
    double q2 = xi * m_from_2;
    double Lambdaxi = Lambda(xi);
    double Gmxi = Gm(xi);
    double xi2 = sqr(xi);
    double xi3 = xi2 * xi;
    switch (channel) {
    case MesonChannel::Scalar : {
        double term1 = sqr(ffp(q2)) * cube(Lambdaxi) / 3. / xi3;
        double term2 = sqr(ffp(q2)) * Lambdaxi * Gmxi * lambda(1., mr_h, xi) / 2. / xi3;
        double term3 = sqr(ff0(q2)) * Lambdaxi * Gmxi * sqr(1. - mr_h) / 2. / xi3;
        return term1 + term2 + term3;
    }
    case MesonChannel::Vector : {
        double fq2 = ff(q2);
        double Fxi = F(xi);
        double Gpxi = Gp(xi);
        double apq2 = ap(q2);
        double amq2 = am(q2);
        double term1 = m_from_2 * mr_h / 3. / xi2 * g(q2) * Lambdaxi * Fxi * (2 * xi2 - Gpxi);
        double term2 = 1. / 24. / m_from_2 / xi3 * sqr(fq2) * Lambdaxi * (3 * Fxi * (xi2 - sqr(mr_l - mr_N)) - sqr(Lambdaxi) + 12 * mr_h * xi * (2 * xi2 - Gpxi));
        double term3 = m_from_2 / 24. / xi3 * sqr(apq2) * Lambdaxi * Fxi * (Fxi * (2 * xi2 - Gpxi) + 3 * Gmxi * sqr(1 - mr_h));
        double term4 = m_from_2 / 8. / xi * sqr(amq2) * Lambdaxi * Fxi * Gmxi;
        double term5 = 1. / 12. / xi3 * fq2 * apq2 * (3 * xi * Fxi * Gmxi + (1 - xi - mr_h) * (3 * Fxi * (xi2 - sqr(mr_l - mr_N)) - sqr(Lambdaxi)));
        double term6 = 1. / 4. / xi2 / fq2 * amq2 * Lambdaxi * Gmxi;
        double term7 = m_from_2 / 4. / xi2 * apq2 * amq2 * Lambdaxi * Fxi * Gmxi * (1 - mr_h);
        return term1 + term2 + term3 + term4 + term5 + term6 + term7;
    }
    default : std::stringstream mess;
//         mess << " unknown decay channel fnSwitch = " << counter;
        info->errorMsg("Warning in StauWidths::function:", mess.str());
        return 0.;
    }
}

}
