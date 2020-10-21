#include "pythia.hh"

#include "ResonanceWidths.hh"
#include "id.hh"
#include "io.hh"
#include "string.hh"

namespace hnl {

namespace {

const bool debug = true;

enum class quark {up, down, strange, charm, bottom, top};

auto& operator<<(std::ostream& stream, quark id) {
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

std::pair<quark, quark> quarks(int id) {
    switch (id) {
        case 211 : return {quark::up, quark::down};
        case 321 : return {quark::up, quark::strange};
        case 411 : return {quark::charm, quark::down};
        case 431 : return {quark::charm, quark::strange};
        case 521 : return {quark::up, quark::bottom};
        case 541 : return {quark::charm, quark::bottom};
        case 111 : return {quark::down, quark::down};
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

bool is_up_type(quark q) {
    return q == quark::up || q == quark::charm || q == quark::top;
}

bool is_down_type(quark q) {
    return q == quark::down || q == quark::strange || q == quark::bottom;
}

struct BRatio {
    BRatio(int a, double b, int c) : onMode(a), bRatio(b), meMode(c) {};
    int onMode;
    double bRatio;
    int meMode;
};

}

bool MesonResonance::can_two_body() {
    if (idRes == 553) return true;
    if (idRes == 443) return true;
    if (idRes == 431) return true;
    auto qs = quarks(idRes);
    return (is_down_type(qs.first) && is_up_type(qs.second)) || (is_up_type(qs.first) && is_down_type(qs.second));
}

bool MesonResonance::can_three_body(int meson) {
    if (idRes == meson) return false;
    if (idRes == 443) return false;
    if (idRes == 411 && meson == 213) return false;
    if (idRes == 431 && meson == 313) return false;
    if (idRes == 431 && meson == 323) return false;
    if (idRes == 431 && meson == 311) return false;
    if (idRes == 431 && meson == 421) return false;
    if (meson == 221) return idRes == 431 ? true : false;
    if (idRes == 511 && meson == 313) return false;
    if (idRes == 521 && meson == 323) return false;
    if (idRes == 421 && meson == 213) return false;
    if (idRes == 521 && meson == 111) return true;
    if (idRes == 321 && meson == 111) return true;
    if (particleDataPtr->chargeType(idRes) == particleDataPtr->chargeType(meson)) return false;
    if (particleDataPtr->m0(idRes) < particleDataPtr->m0(meson)) return false;
    if (idRes == 541) return false;
    auto qs_1 = quarks(idRes);
    auto qs_2 = quarks(meson);
    return qs_1.first == qs_2.first || qs_1.first == qs_2.second || qs_1.second == qs_2.second || qs_1.second == qs_2.first;
}

MesonResonance::MesonResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int id_from) :
    neutrino_coupling(neutrino_coupling_) {
    if (debug) print("Meson Resonance for", id_from);
    initBasic(id_from);
    particlePtr = pythia.particleData.particleDataEntryPtr(idRes);
    std::vector<BRatio> vector;
    for_each(*particlePtr, [&vector](Pythia8::DecayChannel & channel) {
        channel.meMode(101);
        vector.emplace_back(channel.onMode(), channel.bRatio(), 101);
    });
// for (auto pos = 0; pos < particlePtr->sizeChannels(); ++pos) {
// particlePtr->channel(pos).meMode(101);
// vector.emplace_back(particlePtr->channel(pos).onMode(), particlePtr->channel(pos).bRatio(), 101);
// }
    init(&pythia.info, &pythia.settings, &pythia.particleData, pythia.couplingsPtr);
    for (std::size_t pos = 0; pos < vector.size(); ++pos) {
        particlePtr->channel(pos).bRatio(vector.at(pos).bRatio);
        particlePtr->channel(pos).meMode(vector.at(pos).meMode);
        particlePtr->channel(pos).onMode(vector.at(pos).onMode);
    }
}

bool MesonResonance::initBSM() {
    if (debug) print("initBSM");
    if (idRes == 311) print("do not use", 311);
    if (particlePtr->mWidth() == 0. || particlePtr->tau0() > 0.) particlePtr->setMWidth(tau_to_Gamma(particlePtr->tau0()));
    else if (particlePtr->mWidth() > 0. || particlePtr->tau0() == 0.) particlePtr->setTau0(Gamma_to_tau(particlePtr->mWidth()));
    else print("that was unexpectet");
    if (debug) print("Particle", idRes, particlePtr->name(), "mass", particlePtr->m0(), "tau", particlePtr->tau0(), particlePtr->mWidth());
    particlePtr->setIsResonance(false);
    particlePtr->setMayDecay(true);
    return true;
}

void MesonResonance::add_two_body() {
    if (!can_two_body()) return;
    for (auto neutrino : heavy_neutral_leptons()) for (auto lepton : charge_leptons()) particlePtr->addChannel(1, 0., 0, neutrino, -lepton);
}

void MesonResonance::add_three_body(int meson) {
    if (!can_three_body(meson)) return;
    for (auto neutrino : heavy_neutral_leptons()) for (auto lepton : charge_leptons()) particlePtr->addChannel(1, 0., 0, neutrino, particlePtr->chargeType() == 0 ? lepton : - lepton, meson);
}

std::vector<int> MesonResonance::mesons() {
    return {111, 211, 311, 321, 411, 421, 431, 113, 213, 313, 323, 413, 423, 433, 221};
}

bool MesonResonance::allowCalc() {
    if (debug) print("allowCalc");
    particlePtr->clearChannels();
    add_two_body();
    for (auto meson : mesons()) add_three_body(meson);
    return false;
    return true;
}

void MesonResonance::initConstants() {
    if (debug) print("initConstants");
    three_body_width.set_pointers(particleDataPtr);
}

double decay_constant(int id) { // GeV
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
        case 213 : return 0.162; // GeV^2
        case 413 : return 0.535;
        case 433 : return 0.650;
        case 113 : return 0.162;
        case 223 : return 0.153;
        case 333 : return 0.234;
        case 443 : return 1.29;
        default : print("decay_constant for meson", id, "not known");
    }
    return 0.;
}

double NeutrinoResonance::correction_factor(int id) {
    switch (id) {
        case 113 : return 1 - 2 * couplingsPtr->sin2thetaW();
        case 223 : return couplingsPtr->sin2thetaW() * 4 / 3;
        case 333 : return couplingsPtr->sin2thetaW() * 4 / 3 - 1;
        case 443 : return 1 - couplingsPtr->sin2thetaW() * 8 / 3;
        default : print("correction_factor for meson", id, "not known");
    }
    return 0.;
}

std::pair<int, int> quark_pair(int id) {
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
        case 213 : return {1, 1};
        case 413 : return {2, 1};
        case 433 : return {2, 2};
        default : print("quark pair for meson", id, "not known");
    }
    return {0, 0};
}

double clebsch_gordan_2(int id) {
    return id == 111 || id == 113 ? .5 : 1;
}

double MesonResonance::CKM2(int id) {
    std::pair<int, int> pair = quark_pair(id);
    return couplingsPtr->V2CKMgen(pair.first, pair.second);
}

auto get_up_type(std::array<quark, 4> const& quarks) {
    std::map<quark, int> map;
    for (auto quark : quarks) if (is_up_type(quark)) ++map[quark];
    for (auto pair : map) if (pair.second == 1 || pair.second == 3) return pair.first;
    print("up type", "do not end up here", quarks);
    return quark::top;
}

auto get_down_type(std::array<quark, 4> const& quarks) {
    std::map<quark, int> map;
    for (auto quark : quarks) if (is_down_type(quark)) ++map[quark];
    for (auto pair : map) if (pair.second == 1 || pair.second == 3) return pair.first;
    print("no result for", quarks);
    return quark::top;
}

int up_idx(quark up) {
    switch (up) {
        case quark::up : return 1;
        case quark::charm : return 2;
        case quark::top : return 3;
        default : print(up, "is not an up type");
    }
    return 0;
}

int up_idx(int up) {
    switch (up) {
        case 2 : return 1;
        case 4 : return 2;
        case 6 : return 3;
        default : print(up, "is not an up type");
    }
    return 0;
}

int down_idx(quark down) {
    switch (down) {
        case quark::down : return 1;
        case quark::strange : return 2;
        case quark::bottom : return 3;
        default : print(down, "is not an down type");
    }
    return 0;
}

int down_idx(int down) {
    switch (down) {
        case 1 : return 1;
        case 3 : return 2;
        case 5 : return 3;
        default : print(down, "is not an down type");
    }
    return 0;
}

double MesonResonance::CKM2(int id_1, int id_2) {
    if (id_1 == 431 && id_2 == 221) return couplingsPtr->VCKMgen(2, 2);
    auto pair_1 = quarks(id_1);
    auto pair_2 = quarks(id_2);
    auto up = get_up_type({pair_1.first, pair_1.second, pair_2.first, pair_2.second});
    auto down = get_down_type({pair_1.first, pair_1.second, pair_2.first, pair_2.second});

    if (debug) print("In", idRes, id1Abs, id2Abs, id3Abs, "take CKM", up_idx(up), down_idx(down));
    return couplingsPtr->V2CKMgen(up_idx(up), down_idx(down));
}

void MesonResonance::calcPreFac(bool) {
    if (debug) print("calcPreFac");
}

namespace {

double lambda(double a, double b, double c) {
    auto res = sqr(a) + sqr(b) + sqr(c) - 2. * a * b - 2. * a * c - 2. * b * c;
    if (res < 0) print("lambda is neg");
    return res;
}

}

void check_channel(Pythia8::DecayChannel const& channel) {
    if (channel.product(0) != 11 || channel.product(1) != -11) print("should be 11", channel.product(0), channel.product(1));
}

void MesonResonance::calcWidth(bool test) {
    if (debug) print("calcWidth", idRes, id1Abs, id2Abs, id3Abs, "with mass", particleDataPtr->m0(id1Abs));
    print("multiplicity", mult);
    auto id_lep = mult == 2 ? id2Abs : id3Abs;
    preFac = neutrino_coupling(id1Abs, id_lep + 1) * sqr(couplingsPtr->GF()) * Pythia8::pow3(mHat) / 8. / M_PI;
    if (debug) print("preFac", preFac, neutrino_coupling(id1Abs, id_lep + 1), sqr(couplingsPtr->GF()), Pythia8::pow3(mHat));
    widNow = 0.;
    if (preFac <= 0.) {
        print("zero prefactor", idRes, id1Abs, id2Abs, id3Abs);
        particlePtr->setHasChanged(sum > 0.);
        return;
    }
    if (!is_heavy_neutral_lepton(id1Abs)) {
        print(test, "The first particle should be a neutrino", idRes, id1Abs, id2Abs, id3Abs);
        particlePtr->setHasChanged(sum > 0.);
        return;
    }
    if (mHat < mf1 + mf2 + (mult == 2 ? 0 : mf3)) {
        if (debug) print("no phase space", idRes, "to", id1Abs, id2Abs, id3Abs, "mass", mHat,  mf1, mf2, mf3, "sum", mf1 + mf2 + mf3, mf1 + mf2);
        particlePtr->setHasChanged(sum > 0.);
        return;
    }

    if (idRes == 24) { // Not used at the moment
        print("do the W");
        auto thetaWRat = 1. / (12. * couplingsPtr->sin2thetaW());
        alpEM = couplingsPtr->alphaEM(sqr(mHat));
        alpS = couplingsPtr->alphaS(sqr(mHat));
        colQ = 3. * (1. + alpS / M_PI);
        preFac = alpEM * thetaWRat * mHat;
        print(alpEM, alpS, colQ);
        if (ps == 0.) return;
        if (is_heavy_neutral_lepton(id1Abs)) {
            auto id_lep = mult == 2 ? id2Abs : id3Abs;
            auto bL = 1.18921 * particleDataPtr->m0(23) * std::sqrt(couplingsPtr->GF()) * neutrino_coupling(id1Abs, id_lep + 1);
            widNow = mRes * std::sqrt(lambda(1, mf1, mf2)) / 24. / M_PI / mRes * sqr(bL) * ((2. - sqr(mr1 - mr2) - mr1 - mr2) - 6. * mr1 * mr2);
            return;
        }
        if ((id1Abs > 5 && id1Abs < 11) || id1Abs > 16) return;
        widNow = preFac * ps * (1. - 0.5 * (mr1 + mr2) - 0.5 * sqr(mr1 - mr2));
        if (id1Abs < 6) widNow *= colQ * couplingsPtr->V2CKMid(id1Abs, id2Abs);
        return;
    }

    switch (mult) {
        case 2 :
            if (is_heavy_neutral_lepton(id1Abs) && id2Abs > 10 && id2Abs < 17) {
                preFac *= (mr1 + mr2 - sqr(mr2 - mr1)) * std::sqrt(lambda(1., mr1, mr2));
                if (debug) print("preFac", preFac);
                if (idRes == 443 || idRes == 553) {
                    preFac *= 27. * mHat / 8. / M_PI / sqr(couplingsPtr->alphaEM(sqr(mHat)));
                    auto term = 4. / 3. * couplingsPtr->sin2thetaW();
                    if (idRes == 443) {
                        preFac *= sqr(1. - 2 * term) / 4.;
                        auto channel = particlePtr->channel(1);
                        check_channel(channel);
                        widNow = preFac * channel.bRatio() * particlePtr->mWidth() ;
                    } else if (idRes == 553) {
                        preFac *= sqr(1. - term);
                        auto channel = particlePtr->channel(2);
                        check_channel(channel);
                        widNow = preFac * channel.bRatio() * particlePtr->mWidth();
                    }
                } else {
                    preFac *= CKM2(idRes);
                    widNow = preFac * sqr(decay_constant(idRes));
                }
            } else print("Two-body for", id1Abs, "and", id2Abs, "not implemented");
            break;
        case 3 :
            if (is_heavy_neutral_lepton(id1Abs) && id3Abs > 10 && id3Abs < 17 && id2Abs > 20) {
                preFac *= sqr(mHat) / 8. / sqr(M_PI) * clebsch_gordan_2(id2Abs) * CKM2(idRes, id2Abs);
                if (is_vector(id2Abs)) preFac *= sqr(mHat) / sqr(mf2);
                if (debug) print("preFac", preFac);
                widNow = preFac * three_body_width.get_width(idRes, id1Abs, id2Abs, id3Abs);
            } else print("Three-body for", id1Abs, ",", id2Abs, ", and", id3Abs, "not implemented");
            break;
        default : print("multiplicity", mult, "not implemented");
    }
    sum += widNow;
//     particlePtr->setHasChanged(sum > 0.);
    if (debug) print("Calculated", widNow, "for", idRes, "to", id2Abs, "and", id1Abs, "and", id3Abs, "with", preFac, "compare to", tau_to_Gamma(particlePtr->tau0()), particlePtr->mWidth(), sum);
}

Pythia8::ParticleDataEntry & MesonResonance::calculate_widths() {
    if (debug) print("calculate widths");
    initBSM();
    allowCalc();
    initConstants();
    calcPreFac();
    mHat = particlePtr->m0();
    double widTot = 0.;
    for (iChannel = 0; iChannel < particlePtr->sizeChannels(); ++iChannel) {
        mult = particlePtr->channel(iChannel).multiplicity();
        calcWidth();
        print("Current width", widNow, "for", idRes, "to", id2Abs, "and", id1Abs, "and", id3Abs, "<------------");
        widTot += widNow;
        particlePtr->channel(iChannel).bRatio(widNow);
    }
    particlePtr->setMWidth(widTot, false);
    particlePtr->setTau0(Gamma_to_tau(widTot));
    particlePtr->rescaleBR();
    return *particlePtr;
}

NeutrinoResonance::NeutrinoResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& coupling, double mass, int id) :
    neutrino_coupling(coupling) {
    if (debug) print("NeutrinoResonance", id);
    initBasic(id);
    particlePtr = pythia.particleData.particleDataEntryPtr(idRes);
    particlePtr->setM0(mass);
    particlePtr->setMMin(mass / 2);
    particlePtr->setMMax(mass * 2);
    init(&pythia.info, &pythia.settings, &pythia.particleData, pythia.couplingsPtr);
    if (debug) print("Particle", idRes, particlePtr->name(), "mass", particlePtr->m0(), "tau", particlePtr->tau0(), particlePtr->mWidth());
}

bool NeutrinoResonance::initBSM() {
    if (debug) print("initBSM");
    particlePtr->setIsResonance(false);
    particlePtr->setMayDecay(true);
    return true;
}

bool NeutrinoResonance::can_two_body(int meson) {
    if (particlePtr->m0() < particleDataPtr->m0(meson)) return false;
    return true;
}

bool NeutrinoResonance::can_three_body() {
    return true;
}

void NeutrinoResonance::add_three_body() {
    if (!can_three_body()) return;
    for (auto lepton : charge_leptons()) {
        if (mHat >= 1) for (auto up : up_type()) for (auto down : down_type()) particlePtr->addChannel(1, 0., 0, lepton, up, -down);
        for (auto neutrino : neutral_leptons()) for (auto lepton_2 : charge_leptons()) if (lepton_2 != lepton && lepton_2 + 1 == neutrino) particlePtr->addChannel(1, 0., 0, lepton, neutrino, -lepton_2);
    }
    for (auto neutrino : neutral_leptons()) {
        for (auto neutrino_2 : neutral_leptons()) particlePtr->addChannel(1, 0., 0, neutrino, neutrino_2, neutrino_2);
        for (auto lepton : charge_leptons()) particlePtr->addChannel(1, 0., 0, neutrino, lepton, -lepton);
        if (mHat >= 1) for (auto quark : quarks()) particlePtr->addChannel(1, 0., 0, neutrino, quark, -quark);
    }
}

void NeutrinoResonance::add_two_body(int meson) {
    if (!can_two_body(meson)) return;
    if (particleDataPtr->chargeType(meson) == 0) for (auto lepton : neutral_leptons()) particlePtr->addChannel(1, 1., 0, lepton, meson);
    else for (auto lepton : charge_leptons()) particlePtr->addChannel(1, 1., 0, lepton, meson);
}

std::vector<int> NeutrinoResonance::mesons() {
    return {111, 211, 321, 411, 431, 113, 213, 413, 433, 221};
}

bool NeutrinoResonance::allowCalc() {
    if (debug) print("allowCalc");
    particlePtr->clearChannels();
    if (mHat < 1) for (auto meson : mesons()) add_two_body(meson);
    add_three_body();
    return true;
}

void NeutrinoResonance::initConstants() {
    if (debug) print("initConstants");
    three_body_width.set_pointers(particleDataPtr, settingsPtr);
}

double NeutrinoResonance::CKM2(int meson) {
    std::pair<int, int> pair = quark_pair(meson);
    return couplingsPtr->V2CKMgen(pair.first, pair.second);
}

double NeutrinoResonance::CKM2(int up, int down) {
    return couplingsPtr->V2CKMgen(up_idx(up), down_idx(down));
}

double NeutrinoResonance::NW(int up, int down) {
    return is_quark(up) && is_quark(down) ? 3. * CKM2(up, down) : 1.;
}

double L(double x2) {
    if (x2 < 0.000026112) return -21.; // The log becomes complex for electrons; this approximate value is irrelevant as it will be multiplied by a very small number.
    auto sqrt = std::sqrt(1. - 4. * x2);
    return std::log((1. - 3. * x2 - (1. - x2) * sqrt) / x2 / (1. + sqrt));
}

double NeutrinoResonance::Cf1(int neutrino, int fermion) {
    double first = is_quark(fermion) ? (is_up_type(fermion) ? 2. : 1.) / 3. : (fermion + 1 == neutrino ? -1 : 1);
    double second = is_quark(fermion) ? (is_up_type(fermion) ? 4. : 1.) / 9. : 1.;
    return 1. / 4. * (1. - 4. * first * couplingsPtr->sin2thetaW() + 8. * second * sqr(couplingsPtr->sin2thetaW()));
}

double NeutrinoResonance::Cf2(int neutrino, int fermion) {
    double first = is_quark(fermion) ? (is_up_type(fermion) ? 1. : 1. / 2.) / 3. : 1. / 2.;
    double second = is_quark(fermion) ? (is_up_type(fermion) ? 4. : 2.) / 3. : 2.;
    double third = is_quark(fermion) ? 1. : (fermion + 1 == neutrino ? -1. : 1.);
    return first * couplingsPtr->sin2thetaW() * (second * couplingsPtr->sin2thetaW() - third);
}

double NeutrinoResonance::NZ(int fermion) {
    return is_quark(fermion) ? 3 : 1;
}

void NeutrinoResonance::calcPreFac(bool) {
    if (debug) print("calcPreFac");
}

double NeutrinoResonance::get_mass(int id) {
    switch (id) {
        case 1 : return settingsPtr->parm("ParticleData:mdRun");
        case 2 : return settingsPtr->parm("ParticleData:muRun");
        case 3 : return settingsPtr->parm("ParticleData:msRun");
        case 4 : return settingsPtr->parm("ParticleData:mcRun");
        case 5 : return settingsPtr->parm("ParticleData:mbRun");
        case 6 : return settingsPtr->parm("ParticleData:mtRun");
        default : return particleDataPtr->m0(id);
    }
}

void NeutrinoResonance::calcWidth(bool) {
    if (debug) print("calcWidth", idRes, id1Abs, id2Abs, id3Abs);

    auto id1Abs = std::abs(particlePtr->channel(iChannel).product(0));

    auto id2Abs = std::abs(particlePtr->channel(iChannel).product(1));
    auto id3Abs = std::abs(particlePtr->channel(iChannel).product(2));

    auto mf1 = get_mass(id1Abs);
    auto mf2 = get_mass(id2Abs);
    auto mf3 = get_mass(id3Abs);

    auto mr1 = sqr(mf1 / mHat);
    auto mr2 = sqr(mf2 / mHat);
    if (!is_lepton(id1Abs)) {
        print("the second daughter should be a lepton", idRes, "to", id1Abs, id2Abs, mult > 2 ? id3Abs : 0);
        return;
    }
    widNow = 0.;
    preFac = sqr(couplingsPtr->GF()) * Pythia8::pow3(mHat) / 8. / M_PI;
    preFac *= neutrino_coupling(idRes, id1Abs % 2 == 0 ? id1Abs : id1Abs + 1);

    if (preFac <= 0.) {
        if (debug) print("preFac", preFac, idRes, "to", id1Abs, id2Abs, mult > 2 ? id3Abs : 0);
        return;
    }

    if (mHat < mf1 + mf2 + (mult == 2 ? 0 : mf3)) {
        if (debug) print("no phase space", idRes, "to", id1Abs, id2Abs, id3Abs, "mass", mHat, "sum", mf1 + mf2 + mf3, mf1 + mf2);
        return;
    }

    switch (mult) {
        case 2 :
            if (is_meson(id2Abs)) {
                preFac *= sqr(decay_constant(id2Abs)) / 2;
                if (is_charge_lepton(id1Abs)) {
                    preFac *= CKM2(id2Abs);
                    if (is_vector(id2Abs)) {
                        preFac /= sqr(mf2); // dimension of form factor
                        widNow = preFac * (sqr(1 - mr1) + mr2 * (1 + mr1) - 2 * sqr(mr2)) * std::sqrt(lambda(1, mr2, mr1));
                    } else {
                        widNow = preFac * (sqr(1 - mr1) - mr2 * (1 + mr1)) * std::sqrt(lambda(1, mr2, mr1));
                    }
                } else if (is_light_neutrino(id1Abs)) {
                    preFac /= 2;
                    if (is_vector(id2Abs)) {
                        preFac *= sqr(correction_factor(id2Abs)) / sqr(mr2); // dimension of form factor
                        widNow = preFac * (1 + 2 * mr2) * sqr(1 - mr2);
                    } else {
                        widNow = preFac * sqr(1 - mr2);
                    }
                } else print("Two-body not implemented", id1Abs, id2Abs, id3Abs);
            } else print("Two-body not implemented", id1Abs, id2Abs, id3Abs);
            break;
        case 3 :
            preFac *= sqr(mHat) / 24. / sqr(M_PI);
            if (is_charge_lepton(id1Abs)) {
                preFac *= NW(id2Abs, id3Abs);
                if (debug) print("preFac", preFac);
                widNow = preFac * three_body_width.get_width(idRes, id1Abs, id2Abs, id3Abs);
            } else if (is_light_neutrino(id1Abs) && id2Abs == id3Abs) {
                if (is_light_neutrino(id2Abs)) {
                    preFac /= 4.;
                    widNow = preFac * (1. + (id1Abs == id2Abs ? 1. : 0.));
                } else {
                    preFac *= NZ(id2Abs);
// if (mr2 <= 1E-8) print("bad mass", mr2, idRes, id1Abs, id2Abs, id3Abs);
                    auto term_1 = (1. - 14. * mr2 - 2. * sqr(mr2) - 12. * cube(mr2)) * std::sqrt(1. - 4. * mr2) + 12. * sqr(mr2) * (sqr(mr2) - 1.) * L(mr2);
                    auto term_2 = mr2 * (2. + 10. * mr2 - 12 * sqr(mr2)) * std::sqrt(1. - 4. * mr2) + 6 * sqr(mr2) * (1. - 2. * mr2 + 2 * sqr(mr2)) * L(mr2);
                    widNow = preFac * (Cf1(id1Abs, id2Abs) * term_1 + 4. * Cf2(id1Abs, id2Abs) * term_2);
                }
            } else print("this is unexpectet");
            break;
        default : print("multiplicity", mult, "not implemented");
    }
    if (is_quark(id2Abs) && is_quark(id3Abs)) {
        alpS = couplingsPtr->alphaS(sqr(mHat));
        auto mod = alpS / M_PI;
        widNow *= 1. + mod + 5.2 * sqr(mod) + 26.4 * cube(mod);
    }
    if (debug) sum += widNow;
    if (debug) print("Calculated", widNow, "for", idRes, "to", id1Abs, "and", id2Abs, "and", id3Abs, "with", preFac, "compare to", tau_to_Gamma(particlePtr->tau0()), particlePtr->mWidth(), sum);
    if (debug) print(minWidth, particlePtr->isResonance(), particlePtr->mayDecay(), particlePtr->doExternalDecay(), particlePtr->isVisible(), particlePtr->doForceWidth(), particlePtr->mWidth());
    if (debug) print(sum, Gamma_to_tau(sum));
}

Pythia8::ParticleDataEntry & NeutrinoResonance::calculate_widths() {
    if (debug) print("calculate widths");
    initBSM();
    allowCalc();
    initConstants();
    calcPreFac();
    mHat = particlePtr->m0();
    double widTot = 0.;
    for (iChannel = 0; iChannel < particlePtr->sizeChannels(); ++iChannel) {
        mult = particlePtr->channel(iChannel).multiplicity();
        calcWidth();
        widTot += widNow;
        particlePtr->channel(iChannel).bRatio(widNow);
    }
    particlePtr->setMWidth(widTot, false);
    particlePtr->setTau0(Gamma_to_tau(widTot));
    particlePtr->rescaleBR();
    return *particlePtr;
}

}
