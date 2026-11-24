#include <cmath>
#include <complex>
#include <vector>

#include "Pythia8/Pythia.h"

#include "BaryonDecays.hh"
#include "id.hh"
#include "math.hh"
#include "pythia.hh"

namespace hnl {

namespace {

// --- Dirac matrices in the Dirac representation, as 4 x 4 complex arrays

using Complex = std::complex<double>;
using Matrix = std::array<Complex, 16>;

Matrix operator*(Matrix const& a, Matrix const& b) {
    Matrix c{};
    for (int i = 0; i < 4; ++i) for (int k = 0; k < 4; ++k) if (a[4 * i + k] != 0.)
        for (int j = 0; j < 4; ++j) c[4 * i + j] += a[4 * i + k] * b[4 * k + j];
    return c;
}

Matrix operator+(Matrix a, Matrix const& b) {
    for (int i = 0; i < 16; ++i) a[i] += b[i];
    return a;
}

Matrix operator-(Matrix a, Matrix const& b) {
    for (int i = 0; i < 16; ++i) a[i] -= b[i];
    return a;
}

Matrix operator*(Complex x, Matrix a) {
    for (auto& element : a) element *= x;
    return a;
}

Complex trace(Matrix const& a) {
    return a[0] + a[5] + a[10] + a[15];
}

Matrix make(std::array<Complex, 16> const& elements) {
    return elements;
}

const Complex I(0., 1.);
const Matrix one = make({1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1});
const std::array<Matrix, 4> gamma{
    make({1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1}),
    make({0, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, 0, 0}),
    make({0, 0, 0, -I, 0, 0, I, 0, 0, I, 0, 0, -I, 0, 0, 0}),
    make({0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 1, 0, 0})};
const Matrix gamma5 = make({0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0});
const std::array<double, 4> metric{1., -1., -1., -1.};

// gamma^0 A^dagger gamma^0
Matrix bar(Matrix const& a) {
    Matrix dagger;
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) dagger[4 * i + j] = std::conj(a[4 * j + i]);
    return gamma[0] * dagger * gamma[0];
}

double component(Pythia8::Vec4 const& v, int mu) {
    switch (mu) {
        case 0 : return v.e();
        case 1 : return v.px();
        case 2 : return v.py();
        default : return v.pz();
    }
}

Matrix slash(Pythia8::Vec4 const& v) {
    return Complex(v.e()) * gamma[0] - Complex(v.px()) * gamma[1] - Complex(v.py()) * gamma[2] - Complex(v.pz()) * gamma[3];
}

double kallen(double a, double b, double c) {
    return sqr(a) + sqr(b) + sqr(c) - 2. * a * b - 2. * a * c - 2. * b * c;
}

// --- form factors of arXiv:1503.01421: nominal fit in the physical limit, pole masses and t_+ of its table of pole
//     masses, parameters from its ancillary file LambdabLambdac_results.dat

const double m_Bc = 6.276; // GeV, as there

double z_expansion(double q2, double a_0, double a_1, double pole, double m_b, double m_c) {
    double t_plus = sqr(pole); // t_+ = m_pole^2 for Lambda_b -> Lambda_c
    double t_0 = sqr(m_b - m_c);
    double z = (std::sqrt(t_plus - q2) - std::sqrt(t_plus - t_0)) / (std::sqrt(t_plus - q2) + std::sqrt(t_plus - t_0));
    return (a_0 + a_1 * z) / (1. - q2 / sqr(pole));
}

// masses at which the lattice form factors are evaluated (their q2 uses the experimental baryon masses)
const double m_Lambda_b = 5.61960;
const double m_Lambda_c = 2.28646;

// Gauss-Legendre nodes and weights on [0, 1]
std::pair<std::vector<double>, std::vector<double>> gauss_legendre(int n) {
    std::vector<double> x(n), w(n);
    for (int i = 0; i < n; ++i) {
        double z = std::cos(M_PI * (i + 0.75) / (n + 0.5));
        double derivative = 0.;
        for (int iteration = 0; iteration < 100; ++iteration) {
            double p0 = 1., p1 = z;
            for (int k = 2; k <= n; ++k) {
                double p2 = ((2. * k - 1.) * z * p1 - (k - 1.) * p0) / k;
                p0 = p1;
                p1 = p2;
            }
            derivative = n * (z * p1 - p0) / (sqr(z) - 1.);
            double step = p1 / derivative;
            z -= step;
            if (std::abs(step) < 1e-15) break;
        }
        x[i] = (1. - z) / 2.;
        w[i] = 1. / ((1. - sqr(z)) * sqr(derivative));
    }
    return {x, w};
}

// momenta in the Lambda_b rest frame for q2 = (k_l + k_N)^2 and the polar angle of l in the l N frame
void kinematics(double m_b, double m_c, double m_l, double m_N, double q2, double cos_theta, Pythia8::Vec4& p_c, Pythia8::Vec4& k_l, Pythia8::Vec4& k_N) {
    double p_l = std::sqrt(kallen(q2, sqr(m_l), sqr(m_N)) / q2) / 2.;
    double p_q = std::sqrt(kallen(sqr(m_b), q2, sqr(m_c))) / 2. / m_b;
    double sin_theta = std::sqrt(std::max(0., 1. - sqr(cos_theta)));
    k_l = Pythia8::Vec4(p_l * sin_theta, 0., p_l * cos_theta, std::sqrt(sqr(p_l) + sqr(m_l)));
    k_N = Pythia8::Vec4(-k_l.px(), 0., -k_l.pz(), std::sqrt(sqr(p_l) + sqr(m_N)));
    double beta = p_q / std::sqrt(sqr(p_q) + q2);
    k_l.bst(0., 0., beta);
    k_N.bst(0., 0., beta);
    p_c = Pythia8::Vec4(0., 0., -p_q, std::sqrt(sqr(p_q) + sqr(m_c)));
}

}

std::array<double, 6> lambda_b_form_factors(double q2) {
    double pole_f = m_Bc + 0.056, pole_f0 = m_Bc + 0.449, pole_g = m_Bc + 0.492, pole_g0 = m_Bc;
    double a0_gpp = 0.6846557018386; // shared by g_+ and g_perp: g_perp(q2_max) = g_+(q2_max)
    return {z_expansion(q2, 0.81458549275845, -4.8987680713973, pole_f, m_Lambda_b, m_Lambda_c),
            z_expansion(q2, 0.74392566003686, -4.6476511382467, pole_f0, m_Lambda_b, m_Lambda_c),
            z_expansion(q2, 1.077996185173, -6.4170836206134, pole_f, m_Lambda_b, m_Lambda_c),
            z_expansion(q2, a0_gpp, -4.4311552783827, pole_g, m_Lambda_b, m_Lambda_c),
            z_expansion(q2, 0.73960499809758, -4.3664554520299, pole_g0, m_Lambda_b, m_Lambda_c),
            z_expansion(q2, a0_gpp, -4.4633624514401, pole_g, m_Lambda_b, m_Lambda_c)};
}

double lambda_b_matrix_element(Pythia8::Vec4 const& p, Pythia8::Vec4 const& p_c, Pythia8::Vec4 const& k_l, Pythia8::Vec4 const& k_N) {
    double m_b = p.mCalc(), m_c = p_c.mCalc(), m_l = k_l.mCalc(), m_N = k_N.mCalc();
    auto q = p - p_c;
    double q2 = q.m2Calc();
    auto [f_plus, f_0, f_perp, g_plus, g_0, g_perp] = lambda_b_form_factors(q2);
    double s_plus = sqr(m_b + m_c) - q2, s_minus = sqr(m_b - m_c) - q2;
    // <Lambda_c|c gamma^mu (1 - gamma_5) b|Lambda_b> = u_c J^mu u_b, with the helicity form factors as defined in arXiv:1503.01421
    std::array<Matrix, 4> current;
    for (int mu = 0; mu < 4; ++mu) {
        double p_mu = component(p, mu), p_c_mu = component(p_c, mu), q_mu = component(q, mu);
        double longitudinal = p_mu + p_c_mu - (sqr(m_b) - sqr(m_c)) * q_mu / q2;
        auto vector = Complex(f_0 * (m_b - m_c) * q_mu / q2 + f_plus * (m_b + m_c) / s_plus * longitudinal
                              - f_perp * (2. * m_c / s_plus * p_mu + 2. * m_b / s_plus * p_c_mu)) * one + Complex(f_perp) * gamma[mu];
        auto axial = Complex(-1.) * gamma5 * (Complex(g_0 * (m_b + m_c) * q_mu / q2 + g_plus * (m_b - m_c) / s_minus * longitudinal
                     + g_perp * (2. * m_c / s_minus * p_mu - 2. * m_b / s_minus * p_c_mu)) * one + Complex(g_perp) * gamma[mu]);
        current[mu] = vector - axial;
    }
    auto baryon_c = slash(p_c) + Complex(m_c) * one;
    auto baryon_b = slash(p) + Complex(m_b) * one;
    auto lepton = slash(k_l) + Complex(m_l) * one;
    auto neutrino = (one - gamma5) * (slash(k_N) - Complex(m_N) * one) * (one + gamma5);
    std::array<Matrix, 4> left, right;
    for (int mu = 0; mu < 4; ++mu) {
        left[mu] = baryon_c * current[mu] * baryon_b;
        right[mu] = bar(current[mu]);
    }
    Complex sum = 0.;
    for (int mu = 0; mu < 4; ++mu) for (int nu = 0; nu < 4; ++nu) {
        auto hadronic = trace(left[mu] * right[nu]);
        auto leptonic = trace(lepton * gamma[mu] * neutrino * gamma[nu]) * metric[mu] * metric[nu]; // lower indices
        sum += hadronic * leptonic;
    }
    return sum.real() / 4.; // 1/2 from G_F / sqrt(2), 1/2 from the Lambda_b spin
}

double lambda_b_width(double m_b, double m_c, double m_l, double m_N) {
    double q2_min = sqr(m_l + m_N), q2_max = sqr(m_b - m_c);
    if (q2_min >= q2_max) return 0.;
    static auto const nodes_q = gauss_legendre(64);
    static auto const nodes_c = gauss_legendre(24);
    // q2 = q2_min + (q2_max - q2_min) (3 u^2 - 2 u^3), u in [0, 1], smooths the square-root endpoints
    double result = 0.;
    for (std::size_t i = 0; i < nodes_q.first.size(); ++i) {
        double u = nodes_q.first[i];
        double q2 = q2_min + (q2_max - q2_min) * (3. * sqr(u) - 2. * cube(u));
        double jacobian = (q2_max - q2_min) * 6. * u * (1. - u);
        double lambda_q = kallen(sqr(m_b), q2, sqr(m_c)), lambda_l = kallen(q2, sqr(m_l), sqr(m_N));
        if (lambda_q <= 0. || lambda_l <= 0.) continue;
        double inner = 0.;
        for (std::size_t j = 0; j < nodes_c.first.size(); ++j) {
            Pythia8::Vec4 p_c, k_l, k_N;
            kinematics(m_b, m_c, m_l, m_N, q2, 2. * nodes_c.first[j] - 1., p_c, k_l, k_N);
            inner += 2. * nodes_c.second[j] * lambda_b_matrix_element(Pythia8::Vec4(0., 0., 0., m_b), p_c, k_l, k_N);
        }
        result += nodes_q.second[i] * jacobian * std::sqrt(lambda_q * lambda_l) / q2 * inner;
    }
    // dGamma = |M|^2 dm_lN^2 dm_cl^2 / (256 pi^3 m_b^3), dm_cl^2 = sqrt(lambda_q lambda_l) / (2 q2) dcos_theta
    return result / 512. / Pythia8::pow3(M_PI) / Pythia8::pow3(m_b);
}

BaryonResonance::BaryonResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int id_from) :
    neutrino_coupling(neutrino_coupling_) {
    initBasic(id_from);
    particlePtr = pythia.particleData.particleDataEntryPtr(idRes);
    tau0_sm = particlePtr->tau0();
}

bool BaryonResonance::init(Pythia8::Info* infoPtrIn) {
    sum = 0.;
    if (!ResonanceWidths::init(infoPtrIn)) return false;
    // as MesonResonance: width SM plus HNL channels, the channels open for the particle add up to one
    double width = tau_to_Gamma(tau0_sm) + sum;
    particlePtr->rescaleBR(particlePtr->mWidth() / width);
    particlePtr->setMWidth(width, false);
    particlePtr->setTau0(Gamma_to_tau(width));
    return true;
}

bool BaryonResonance::initBSM() {
    particlePtr->setTau0(tau0_sm); // the previous init() has changed it
    particlePtr->setMWidth(tau_to_Gamma(tau0_sm));
    particlePtr->setIsResonance(false);
    particlePtr->setMayDecay(true);
    return true;
}

bool BaryonResonance::allowCalc() {
    if (channels_added) return true; // init() runs in the constructor of Pythia's resonance list and in pythia.init()
    if (idRes == 5122) for (auto neutrino : heavy_neutral_leptons()) for (auto lepton : charge_leptons())
        particlePtr->addChannel(1, 0., 0, neutrino, lepton, 4122);
    channels_added = true;
    return true;
}

void BaryonResonance::calcWidth(bool) {
    auto const& channel = particlePtr->channel(iChannel);
    if (!has_neutrino(channel)) { // SM channel, keeps its width
        widNow = GammaRes * channel.bRatio();
        return;
    }
    int neutrino = std::abs(channel.product(0)), lepton = std::abs(channel.product(1)), baryon = std::abs(channel.product(2));
    double coupling = neutrino_coupling(neutrino, lepton + 1);
    widNow = 0.;
    if (coupling <= 0.) return;
    widNow = sqr(coupSMPtr->GF()) * coupSMPtr->V2CKMgen(2, 3) * coupling
        * lambda_b_width(particleDataPtr->m0(idRes), particleDataPtr->m0(baryon), particleDataPtr->m0(lepton), particleDataPtr->m0(neutrino));
    sum += widNow;
}

}
