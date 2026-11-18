#include "Pythia8/Pythia.h"

#include "MesonDecays.hh"
#include "id.hh"
#include "math.hh"
#include "pythia.hh"

namespace hnl {

namespace {

const bool debug = false;

double kallen(double a, double b, double c) {
    return sqr(a) + sqr(b) + sqr(c) - 2. * a * b - 2. * a * c - 2. * b * c;
}

// epsilon_{mu nu rho sigma} a^mu b^nu c^rho d^sigma with epsilon_{0123} = 1
double epsilon(Pythia8::Vec4 const& a, Pythia8::Vec4 const& b, Pythia8::Vec4 const& c, Pythia8::Vec4 const& d) {
    double m[4][4] = {{a.e(), a.px(), a.py(), a.pz()}, {b.e(), b.px(), b.py(), b.pz()}, {c.e(), c.px(), c.py(), c.pz()}, {d.e(), d.px(), d.py(), d.pz()}};
    double det = 0.;
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) if (j != i) for (int k = 0; k < 4; ++k) if (k != i && k != j) {
        int l = 6 - i - j - k;
        int sign = ((i < j) + (i < k) + (i < l) + (j < k) + (j < l) + (k < l)) % 2 == 0 ? 1 : -1; // parity of the permutation
        det += sign * m[0][i] * m[1][j] * m[2][k] * m[3][l];
    }
    return det;
}

// v with v * a = epsilon(a, b, c, d) for all a
Pythia8::Vec4 epsilon(Pythia8::Vec4 const& b, Pythia8::Vec4 const& c, Pythia8::Vec4 const& d) {
    double e = epsilon(Pythia8::Vec4(0., 0., 0., 1.), b, c, d);
    double x = epsilon(Pythia8::Vec4(1., 0., 0., 0.), b, c, d);
    double y = epsilon(Pythia8::Vec4(0., 1., 0., 0.), b, c, d);
    double z = epsilon(Pythia8::Vec4(0., 0., 1., 0.), b, c, d);
    return Pythia8::Vec4(-x, -y, -z, e);
}

}

MesonDecays::MesonDecays(Pythia8::Pythia& pythia, std::vector<int> const& mesons_) :
    particle_data(&pythia.particleData), rndm(&pythia.rndm), mesons(mesons_) {
    three_body_width.set_pointers(particle_data);
}

std::vector<int> MesonDecays::handledParticles() {
    return mesons;
}

void MesonDecays::restore() {
    for (auto& [entry, pos, on_mode] : switched_off) entry->channel(pos).onMode(on_mode);
    switched_off.clear();
}

bool MesonDecays::decay(std::vector<int>& ids, std::vector<double>& masses, std::vector<Pythia8::Vec4>& momenta, int line, Pythia8::Event const&) {
    restore();
    int id = ids[0];
    auto entry = particle_data->particleDataEntryPtr(id);
    if (line == forced) for (int pos = 0; pos < entry->sizeChannels(); ++pos) if (!has_neutrino(entry->channel(pos))) {
        switched_off.emplace_back(entry, pos, entry->channel(pos).onMode());
        entry->channel(pos).onMode(0);
    }
    if (!entry->preparePick(id, masses[0])) return false;
    auto const& channel = entry->pickChannel();
    if (!has_neutrino(channel)) { // SM decay: Pythia picks again, without the HNL channels
        for (int pos = 0; pos < entry->sizeChannels(); ++pos) if (has_neutrino(entry->channel(pos))) {
            switched_off.emplace_back(entry, pos, entry->channel(pos).onMode());
            entry->channel(pos).onMode(0);
        }
        return false;
    }
    std::vector<int> products;
    for (int pos = 0; pos < channel.multiplicity(); ++pos) {
        int product = channel.product(pos);
        products.push_back(id < 0 && particle_data->hasAnti(product) ? -product : product);
    }
    int id_N = 0, id_l = 0, id_to = 0;
    for (auto product : products) {
        if (is_heavy_neutral_lepton(std::abs(product))) id_N = product;
        else if (is_lepton(std::abs(product))) id_l = product;
        else id_to = product;
    }
    double m_N = particle_data->m0(id_N);
    double m_l = particle_data->m0(id_l);
    Pythia8::Vec4 k_N, k_l, p_to;
    double m_to = 0.;
    if (products.size() == 2) {
        if (masses[0] < m_N + m_l) return false;
        std::tie(k_N, k_l) = rndm->phaseSpace2(masses[0], m_N, m_l);
    } else {
        bool success = false;
        for (int tries = 0; tries < 10 && !success; ++tries) {
            m_to = particle_data->mSel(id_to);
            success = three_body(id, id_to, masses[0], m_to, m_l, m_N, p_to, k_l, k_N);
        }
        if (!success) return false;
    }
    for (auto product : products) {
        auto momentum = product == id_N ? k_N : product == id_l ? k_l : p_to;
        momentum.bst(momenta[0], masses[0]);
        ids.push_back(product);
        masses.push_back(product == id_N ? m_N : product == id_l ? m_l : m_to);
        momenta.push_back(momentum);
    }
    if (debug) print("decayed", id, "to", products);
    return true;
}

double MesonDecays::matrix_element(int id_from, int id_to, Pythia8::Vec4 const& p, Pythia8::Vec4 const& p_to, Pythia8::Vec4 const& k_l, Pythia8::Vec4 const& k_N) {
    auto P = p + p_to;
    auto q = p - p_to;
    auto form_factors = three_body_width.form_factors(id_from, id_to, q.m2Calc());
    if (!is_vector(std::abs(id_to))) { // eq. 22 of 1805.08567
        double f_minus = (form_factors[1] - form_factors[0]) * (p.m2Calc() - p_to.m2Calc()) / q.m2Calc();
        auto h = form_factors[0] * P + f_minus * q;
        return 2. * (k_l * h) * (k_N * h) - (h * h) * (k_l * k_N);
    }
    // eqs. 23 and 24 of 1805.08567, summed over three real polarisation vectors. The sign of the epsilon term is fixed by
    // V-A: for m_N -> 0 the lepton is harder than the neutrino, as in b -> c l nu where |M|^2 ~ (p_b p_nu)(p_c p_l)
    double result = 0.;
    for (auto e : {Pythia8::Vec4(1., 0., 0., 0.), Pythia8::Vec4(0., 1., 0., 0.), Pythia8::Vec4(0., 0., 1., 0.)}) {
        e.bst(p_to, p_to.mCalc());
        auto a = - form_factors[1] * e - form_factors[2] * (e * p) * P - form_factors[3] * (e * p) * q; // axial current, real
        auto b = form_factors[0] * epsilon(e, P, q); // vector current, imaginary
        result += 2. * ((k_l * a) * (k_N * a) + (k_l * b) * (k_N * b)) - (a * a + b * b) * (k_l * k_N) + 2. * epsilon(a, b, k_l, k_N);
    }
    return result;
}

double MesonDecays::weight(int id_from, int id_to, double m_from, double m_to, double m_l, double m_N, double q2, double cos_theta, double phi, Pythia8::Vec4& p_to, Pythia8::Vec4& k_l, Pythia8::Vec4& k_N) {
    double lambda_q = kallen(sqr(m_from), q2, sqr(m_to));
    double lambda_l = kallen(q2, sqr(m_l), sqr(m_N));
    if (lambda_q <= 0. || lambda_l <= 0.) return 0.;
    double p_l = std::sqrt(lambda_l / q2) / 2.; // in the l N frame
    double p_q = std::sqrt(lambda_q) / 2. / m_from; // in the meson frame
    double sin_theta = std::sqrt(1. - sqr(cos_theta));
    k_l = Pythia8::Vec4(p_l * sin_theta * std::cos(phi), p_l * sin_theta * std::sin(phi), p_l * cos_theta, std::sqrt(sqr(p_l) + sqr(m_l)));
    k_N = Pythia8::Vec4(-k_l.px(), -k_l.py(), -k_l.pz(), std::sqrt(sqr(p_l) + sqr(m_N)));
    double beta = p_q / std::sqrt(sqr(p_q) + q2);
    k_l.bst(0., 0., beta);
    k_N.bst(0., 0., beta);
    p_to = Pythia8::Vec4(0., 0., -p_q, std::sqrt(sqr(p_q) + sqr(m_to)));
    // flat in (m_lN^2, m_hl^2) is flat in (q2, cos_theta) times sqrt(lambda_q lambda_l) / q2
    return matrix_element(id_from, id_to, Pythia8::Vec4(0., 0., 0., m_from), p_to, k_l, k_N) * std::sqrt(lambda_q * lambda_l) / q2;
}

bool MesonDecays::three_body(int id_from, int id_to, double m_from, double m_to, double m_l, double m_N, Pythia8::Vec4& p_to, Pythia8::Vec4& k_l, Pythia8::Vec4& k_N) {
    double q2_min = sqr(m_l + m_N);
    double q2_max = sqr(m_from - m_to);
    if (q2_min >= q2_max) return false;
    auto random_weight = [&]() {
        return weight(id_from, id_to, m_from, m_to, m_l, m_N, q2_min + (q2_max - q2_min) * rndm->flat(), 2. * rndm->flat() - 1., 2. * M_PI * rndm->flat(), p_to, k_l, k_N);
    };
    auto key = std::make_tuple(std::abs(id_from), std::abs(id_to), int(1e4 * m_l), int(1e4 * m_N));
    if (!max_weight.count(key)) {
        double max = 0.;
        for (int point = 0; point < 10000; ++point) max = std::max(max, random_weight());
        max_weight[key] = 1.5 * max;
    }
    auto& max = max_weight[key];
    for (int tries = 0; tries < 100000; ++tries) {
        double now = random_weight();
        if (now > max) {
            print("MesonDecays: weight", now, "above maximum", max, "for", id_from, id_to);
            max = 1.2 * now;
        }
        if (now < max * rndm->flat()) continue;
        double theta = std::acos(2. * rndm->flat() - 1.);
        double phi = 2. * M_PI * rndm->flat();
        for (auto* momentum : {&p_to, &k_l, &k_N}) momentum->rot(theta, phi);
        return true;
    }
    return false;
}

}
