#include "ResonanceWidths.hh"
#include <boost/range/adaptor/indexed.hpp>

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


template<typename Element, template <typename, typename = std::allocator<Element>> class Container>
auto & operator<<(std::ostream& stream, Container<Element> const& container) noexcept
{
    for (auto const& element : boost::adaptors::index(container)) stream << '\n' << element.index() << ": " << element.value();
    return stream;
}

template<typename Element, size_t size>
auto& operator<<(std::ostream& stream, std::array<Element, size> const& container) noexcept
{
    for (auto const& element : boost::adaptors::index(container)) stream << '\n' << element.index() << ": " << element.value();
    return stream;
}

template<typename Object, typename ... Arguments>
void print(Object const& object, Arguments ... arguments) noexcept
{
    std::cout << std::boolalpha << std::scientific << object << ' ';
    print(arguments ...);
}

bool debug = true;

bool is_vector(int id)
{
    return id == 113 || id == 213 || id == 313 || id == 323 || id == 413 || id == 423 || id == 433;
}

}

enum class quark {up, down, strange, charm, bottom, top};

auto& operator<<(std::ostream& stream, quark id) noexcept
{
    switch (id) {
    case quark::up : return stream << "u";
    case quark::down : return stream << "d";
    case quark::strange : return stream << "s";
    case quark::charm : return stream << "c";
    case quark::bottom : return stream << "b";
    case quark::top : return stream << "t";
    default : return stream << "none";
    }
}

std::pair<quark, quark> quarks(int id)
{
    switch (id) {
    case 211 : return {quark::up, quark::down};
    case 321 : return {quark::up, quark::strange};
    case 411 : return {quark::charm, quark::down};
    case 431 : return {quark::charm, quark::strange};
    case 521 : return {quark::up, quark::bottom};
    case 541 : return {quark::charm, quark::bottom};
    case 111 : return {quark::up, quark::down};
    case 221 : return {quark::up, quark::strange};
    case 331 : return {quark::up, quark::strange};
    case 441 : return {quark::charm, quark::strange};

    case 311 : return {quark::down, quark::strange};
    case 511 : return {quark::down, quark::bottom};
    case 421 : return {quark::charm, quark::up};
    case 113 : return {quark::up, quark::up};
    case 531 : return {quark::bottom, quark::strange};
    case 213 : return {quark::up, quark::down};
    case 130 : return {quark::down, quark::strange};
    case 310 : return {quark::down, quark::strange};
    case 553 : return {quark::bottom, quark::bottom};
    case 313 : return {quark::strange, quark::down};
    case 323 : return {quark::strange, quark::up};
    case 413 : return {quark::charm, quark::down};
    case 423 : return {quark::charm, quark::up};
    case 433 : return {quark::charm, quark::strange};
    case 443 : return {quark::charm, quark::charm};
    default : print("meson ", id, " not known", "quarks");
    }
    return {quark::top, quark::top};
}

bool is_up_type(quark q)
{
    return q == quark::up || q == quark::charm || q == quark::top;
}

bool is_down_type(quark q)
{
    return q == quark::down || q == quark::strange || q == quark::bottom;
}

bool MesonResonance::can_two_body()
{
    if (idRes == 553) return true;
    if (idRes == 443) return true;
    auto qs = quarks(idRes);
    return (is_down_type(qs.first) && is_up_type(qs.second)) || (is_up_type(qs.first) && is_down_type(qs.second));
}

bool MesonResonance::can_three_body(int id)
{
    if (idRes == 443) return false;
    if (idRes == id) return false;
    if (idRes == 411 && id == 213) return false;
    if (idRes == 431 && id == 313) return false;
    if (idRes == 431 && id == 323) return false;
    if (idRes == 511 && id == 313) return false;
    if (idRes == 521 && id == 323) return false;
    if (particleDataPtr->chargeType(idRes) == particleDataPtr->chargeType(id)) return false;
    if (particleDataPtr->m0(idRes) < particleDataPtr->m0(id)) return false;
    if (idRes == 541) return false;
    auto qs_1 = quarks(idRes);
    auto qs_2 = quarks(id);
    return qs_1.first == qs_2.first || qs_1.first == qs_2.second || qs_1.second == qs_2.second || qs_1.second == qs_2.first;
}


double tau_to_Gamma(double tau)//mm/c->GeV
{
    return 1.97327E-13 / tau;
}


MesonResonance::MesonResonance(Pythia8::Pythia& pythia, double neutrino_coupling_, int id_from_) :
    neutrino_coupling(neutrino_coupling_)
{
    initBasic(id_from_);
    if (idRes == 311) print("do not use", 311);
    standard_model.init(pythia.settings, &pythia.rndm);
    pythia.particleData.findParticle(idRes)->setMWidth(tau_to_Gamma(pythia.particleData.findParticle(idRes)->tau0()));

//     for (auto pos = 0; pos < particle_data_entry.sizeChannels(); ++pos) {
//         auto channel = particle_data_entry.channel(pos);
//         channels[ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}] = {channel.onMode(), channel.bRatio(), channel.meMode()};
//     }

    for (auto pos = 0; pos < pythia.particleData.findParticle(idRes)->sizeChannels(); ++pos) pythia.particleData.findParticle(idRes)->channel(pos).meMode(101);
    if (debug) print("Particle", idRes, pythia.particleData.findParticle(idRes)->name(), "mass", pythia.particleData.findParticle(idRes)->m0(), "tau", pythia.particleData.findParticle(idRes)->tau0(), pythia.particleData.findParticle(idRes)->mWidth());
}

bool MesonResonance::initBSM()
{
    if (debug) print("initBSM");
    return true;
}

void MesonResonance::AddMissingChannels(Pythia8::ParticleData& particle_data)
{
    auto* particle = particle_data.findParticle(idRes);
//     for (auto const& row : channels) {
//         if (std::get<0>(row.first) == 9900014) continue;
//         if (std::get<1>(row.first) == 9900014) continue;
//         if (std::get<2>(row.first) == 9900014) continue;
//         if (std::get<3>(row.first) == 9900014) continue;
//         if (std::get<4>(row.first) == 9900014) continue;
// //         if (debug)
//             print("add", row.second.onMode, row.second.bRatio, 110, std::get<0>(row.first), std::get<1>(row.first), std::get<2>(row.first), std::get<3>(row.first), std::get<4>(row.first));
//         particle->addChannel(row.second.onMode, row.second.bRatio, row.second.meMode, std::get<0>(row.first), std::get<1>(row.first), std::get<2>(row.first), std::get<3>(row.first), std::get<4>(row.first));
//     }
//     particle->setTau0(particle_data_entry.tau0());
    particle->rescaleBR();
}

bool MesonResonance::allowCalc()
{
    if (debug) print("allowCalc");
    bool done = getChannels();
    std::stringstream idStream;
    idStream << "ID = " << idRes ;
    if (!done) infoPtr->errorMsg("Error in SusyResonanceWidths::allowcalc: " "unable to reset decay table.", idStream.str(), true);
    return done;
}

void MesonResonance::add_two_body(Pythia8::ParticleDataEntry& particle, int neutrino)
{
    if (!can_two_body()) return;
//     particle.chargeType() == 0 ?
//     particle.addChannel(1, 0., 0, 9900014, neutrino) :
    particle.addChannel(1, 0., 0, 9900014, 1 - neutrino);
}

void MesonResonance::add_three_body(Pythia8::ParticleDataEntry& particle, int neutrino, int id)
{
    if (!can_three_body(id)) return;
//     particle.chargeType() == particleDataPtr->chargeType(id) ?
//     particle.addChannel(1, 0., 0, 9900014, neutrino, id) :
    particle.addChannel(1, 0., 0, 9900014, neutrino - 1, id);
}

bool MesonResonance::getChannels()
{
    if (debug) print("getChannels");
    auto& particle = *particleDataPtr->particleDataEntryPtr(idRes);
//     particle.clearChannels();
    std::vector<int> neutrinos{14};
    for (auto const& neutrino : neutrinos) {
        add_two_body(particle, neutrino);
        add_three_body(particle, neutrino, 111);
        add_three_body(particle, neutrino, 211);
        add_three_body(particle, neutrino, 311);
        add_three_body(particle, neutrino, 321);
        add_three_body(particle, neutrino, 411);
        add_three_body(particle, neutrino, 421);
        add_three_body(particle, neutrino, 431);
        add_three_body(particle, neutrino, 113);
        add_three_body(particle, neutrino, 213);
        add_three_body(particle, neutrino, 313);
        add_three_body(particle, neutrino, 323);
        add_three_body(particle, neutrino, 413);
        add_three_body(particle, neutrino, 423);
        add_three_body(particle, neutrino, 433);
    }
    return true;
}

void MesonResonance::initConstants()
{
    if (debug) print("initConstants");
    three_body_width.set_pointers(particleDataPtr);
}

double decay_constant(int id) // GeV
{
    switch (id) {
    case 211 : return 0.1302;
    case 321 : return 0.1556;
    case 411 : return 0.212;
    case 431 : return 0.249;
    case 521 : return 0.187;
    case 541 : return 0.434;
    case 111 : return 0.1302;
    case 221 : return 0.0817;
    case 331 : return -0.0947;
    case 441 : return 0.237;
    default : print("meson ", id, " not known", "decay_constant");
        return 0.;
    }
}

std::pair<int, int> quark_pair(int id)
{
    switch (id) {
    case 211 : return {1, 1};
    case 321 : return {1, 2};
    case 411 : return {2, 1};
    case 431 : return {2, 2};
    case 521 : return {1, 3};
    case 541 : return {2, 3};
    case 111 : return {1, 1};
    case 221 : return {1, 2};
    case 331 : return {1, 2};
    case 441 : return {2, 2};
    default : print("meson ", id, " not known", "quark_pair");
        return {0, 0};
    }
}

double MesonResonance::NeutU(int id_heavy, int id_light)
{
    auto diff = id_heavy - id_light - 9900000;
    if (debug) print("diff", diff);
    double cpl = (diff == 0 || diff == 1) ? neutrino_coupling : 0.;
    if (debug) print("coupl", cpl);
    return cpl;
}

double clepsch_gordan_2(int id)
{
    return id == 111 || id == 113 ? .5 : 1;
}

double MesonResonance::CKM2(int id)
{
    std::pair<int, int> pair = quark_pair(id);
    return standard_model.V2CKMgen(pair.first, pair.second);
}


auto get_up_type(std::array<quark, 4> const& quarks)
{
    quark res = quark::top;
    for (auto quark_1 : quarks) {
        if (is_down_type(quark_1)) continue;
        res = quark_1;
        for (auto quark_2 : quarks) {
            if (is_down_type(quark_2)) continue;
            if (quark_1 != quark_2) return quark_2;
        }
    }
    return res;
}

auto get_down_type(std::array<quark, 4> const& quarks)
{
    quark res = quark::top;
    for (auto quark_1 : quarks) {
        if (is_up_type(quark_1)) continue;
        res = quark_1;
        for (auto quark_2 : quarks) {
            if (is_up_type(quark_2)) continue;
            if (quark_1 != quark_2) return quark_2;
        }
    }
    return res;
}

int up_idx(quark up)
{
    switch (up) {
    case quark::up : return 1;
    case quark::charm : return 2;
    case quark::top : return 3;
    default : print(up, "is not an up type");
    }
    return 0;
}

int down_idx(quark down)
{
    switch (down) {
    case quark::down : return 1;
    case quark::strange : return 2;
    case quark::bottom : return 3;
    default : print(down, "is not an down type");
    }
    return 0;
}

double MesonResonance::CKM2(int id_1, int id_2)
{
    auto pair_1 = quarks(id_1);
    auto pair_2 = quarks(id_2);
    auto up = get_up_type({pair_1.first, pair_1.second, pair_2.first, pair_2.second});
    auto down = get_down_type({pair_1.first, pair_1.second, pair_2.first, pair_2.second});
    return standard_model.V2CKMgen(up_idx(up), down_idx(down));
}

void MesonResonance::calcPreFac(bool) // Common coupling factors.
{
    if (debug) print("calcPreFac");
}

double lambda(double a, double b, double c)
{
    return sqr(a) + sqr(b) + sqr(c) - 2 * a * b - 2 * a * c - 2 * b * c;
}

void MesonResonance::calcWidth(bool)
{
    if (debug) print("calcWidth", idRes, id1Abs, id2Abs, id3Abs);
    if (id1Abs != 9900014 && id2Abs != 9900014 && id3Abs != 9900014) {
        print("do not end up here");
        return;
    }
    if (ps == 0.) {
        if (debug) print("no phase space");
        preFac = 0.;
        widNow = 0.;
        return;
    }
    auto id_lep = mult == 2 ? id2Abs : id3Abs;
    preFac = NeutU(id1Abs, id_lep) * sqr(standard_model.GF()) * Pythia8::pow3(mHat) / 8 / M_PI;
    if (debug) print("preFac", preFac);
    switch (mult) {
    case 2 : if (id1Abs == 9900014 && id2Abs > 10 && id2Abs < 17) {
//             print("ids", idRes, id1Abs, id2Abs);
//             print("masses", mHat, mf1, mf2);
//             print("ratios 1", mr1, mr2);
//             print("ratios 2", sqr(mf1 / mHat), sqr(mf2 / mHat));
//             print("lambda", lambda(1, mr1, mr2));
            preFac *= (mr1 + mr2 - sqr(mr2 - mr1)) * std::sqrt(lambda(1, mr1, mr2));
            if (debug) print("preFac", preFac);
            if (idRes == 443) {
                auto channel = particle_data_entry.channel(1);
                if (channel.product(0) != 11 || channel.product(1) != -11) print("should be 11", channel.product(0), channel.product(1));
                preFac *=  27 * mHat / 32 / M_PI / sqr(standard_model.alphaEM(mHat)) * sqr(1 - 8. / 3. * standard_model.sin2thetaW());
                if (debug) print("preFac", preFac);
                widNow = preFac * channel.bRatio() * (particle_data_entry.mWidth() > 0. ? particle_data_entry.mWidth() > 0. : tau_to_Gamma(particle_data_entry.tau0()));
            } else if (idRes == 553) {
                auto channel = particle_data_entry.channel(2);
                if (channel.product(0) != 11 || channel.product(1) != -11) print("should be 11", channel.product(0), channel.product(1));
                preFac *= 27 * mHat / 8 / M_PI / sqr(standard_model.alphaEM(mHat)) * sqr(-1 + 4. / 3. * standard_model.sin2thetaW());
                if (debug) print("preFac", preFac);
                widNow = preFac * channel.bRatio() * (particle_data_entry.mWidth() > 0. ? particle_data_entry.mWidth() > 0. : tau_to_Gamma(particle_data_entry.tau0()));
            } else {
                preFac *= CKM2(idRes);
                widNow = preFac * sqr(decay_constant(idRes));
            }
        } else print("Two-body not implemented");
        break;
    case 3 : if (id1Abs == 9900014 && id3Abs > 10 && id3Abs < 17 && id2Abs > 20) {
            preFac *= sqr(mHat) / 8 / sqr(M_PI) * clepsch_gordan_2(id2Abs) * CKM2(idRes, id2Abs);
            if (debug) print("preFac", preFac);
            if (is_vector(id2Abs)) preFac *= sqr(mHat) / sqr(mf2);
            widNow = preFac * three_body_width.get_width(idRes, id2Abs, id3Abs);
        } else print("Three-body for", id1Abs, id2Abs, id3Abs, "not implemented");
        break;
    default : print("multiplicity", mult, "not implemented");
    }
    if (debug) print("Calculated", widNow, "for", idRes, "to", id2Abs, "to", id2Abs, "with", preFac, "compare to", tau_to_Gamma(particle_data_entry.tau0()), particle_data_entry.mWidth());
}

}
