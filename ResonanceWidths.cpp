#include "pythia.hh"

#include "ResonanceWidths.hh"
#include "id.hh"
#include "io.hh"
#include "string.hh"

namespace hnl {

namespace {

const bool debug = false;

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
        case 441 : return {quark::charm, quark::charm};

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

}

bool MesonResonance::can_two_body() {
    if (idRes == 553) return true;
    if (idRes == 443) return true;
    if (idRes == 431) return true;
    auto qs = quarks(idRes);
    return (is_down_type(qs.first) && is_up_type(qs.second)) || (is_up_type(qs.first) && is_down_type(qs.second));
}

bool MesonResonance::can_three_body(int meson) {
    if (meson == 443 || meson == 441) return idRes == 541;
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
    tau0_sm = particlePtr->tau0();
    width_sm = particlePtr->mWidth();
}

bool MesonResonance::init(Pythia8::Info* infoPtrIn) {
    sum = 0.;
    if (!ResonanceWidths::init(infoPtrIn)) return false;
    // Pythia lists some decays separately for particle and antiparticle (onMode 2 and 3) and init() adds up all channels.
    // Width: SM plus HNL channels. Branching ratios: the channels open for the particle add up to one, as in Pythia.
    double width = (tau0_sm > 0. ? tau_to_Gamma(tau0_sm) : width_sm) + sum;
    particlePtr->rescaleBR(particlePtr->mWidth() / width);
    particlePtr->setMWidth(width, false);
    particlePtr->setTau0(Gamma_to_tau(width));
    return true;
}

bool MesonResonance::initBSM() {
    if (debug) print("initBSM");
    if (idRes == 311) print("do not use", 311);
    particlePtr->setTau0(tau0_sm); // the previous init() has changed them
    particlePtr->setMWidth(width_sm);
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
    if (particleDataPtr->hasAnti(meson)) meson = -meson; // positive mesons decay via W+, e.g. B+ -> anti-D0 l+ N
    for (auto neutrino : heavy_neutral_leptons()) for (auto lepton : charge_leptons()) particlePtr->addChannel(1, 0., 0, neutrino, - lepton, meson);
}

std::vector<int> MesonResonance::mesons() {
    return {111, 211, 311, 321, 411, 421, 431, 113, 213, 313, 323, 413, 423, 433, 221, 441, 443};
}

bool MesonResonance::allowCalc() {
    if (debug) print("allowCalc");
//     particlePtr->clearChannels();
    if (channels_added) return true; // init() runs twice: in the constructor and in pythia.init()
    add_two_body();
    for (auto meson : mesons()) add_three_body(meson);
    channels_added = true;
    return true;
    return false;
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
        // not in 1805.08567; from tau decays with the l V formula (eq. 40) and Pythia's CKM, as g_rho there:
        case 323 : return 0.1795; // BR(tau -> K* nu) / BR(tau -> pi nu) = 1.20 / 10.82 %
        case 20213 : return 0.2861; // BR(tau -> 3 pi nu) / BR(tau -> pi nu) = (9.31 + 9.26) / 10.82 %, m_a1 = 1.23 GeV
        case 20113 : return 0.2861;
        default : print("decay_constant for meson", id, "not known");
    }
    return 0.;
}

// kappa of <0|J^Z|V> = kappa g_V / sqrt(2) eps (1805.08567): sqrt(2) sum_q c_q (T3_q - 2 Q_q sin^2 theta_W), with the quark
// content c_q of rho0, omega = (u ubar -+ d dbar) / sqrt(2), phi = s sbar, J/psi = c cbar, for which the g_V of 1805.08567
// are the standard decay constants (f m). This gives its kappa for rho0; its table has 2x (omega) and sqrt(2)x (phi, J/psi)
// these values. J/psi -> nu nubar in MesonResonance uses the same (1 - 8/3 sin^2 theta_W) / 2.
double NeutrinoResonance::correction_factor(int id) {
    switch (id) {
        case 113 : return 1 - 2 * coupSMPtr->sin2thetaW();
        case 223 : return coupSMPtr->sin2thetaW() * 2 / 3;
        case 333 : return (1 - coupSMPtr->sin2thetaW() * 4 / 3) / M_SQRT2;
        case 443 : return (1 - coupSMPtr->sin2thetaW() * 8 / 3) / M_SQRT2;
        case 20113 : return 1; // isovector axial current without electromagnetic part (appendix C of 1805.08567)
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
    return coupSMPtr->V2CKMgen(pair.first, pair.second);
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
    if (id_1 == 431 && id_2 == 221) return coupSMPtr->VCKMgen(2, 2);
    auto pair_1 = quarks(id_1);
    auto pair_2 = quarks(id_2);
    auto up = get_up_type({pair_1.first, pair_1.second, pair_2.first, pair_2.second});
    auto down = get_down_type({pair_1.first, pair_1.second, pair_2.first, pair_2.second});

    if (debug) print("In", idRes, id1Abs, id2Abs, id3Abs, "take CKM", up_idx(up), down_idx(down));
    return coupSMPtr->V2CKMgen(up_idx(up), down_idx(down));
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

// branching ratio of the SM channel V -> e+ e-
double electron_branching_ratio(Pythia8::ParticleDataEntry& particle) {
    double result = 0.;
    for_each(particle, [&result](Pythia8::DecayChannel const& channel) {
        if (channel.multiplicity() == 2 && channel.product(0) == 11 && channel.product(1) == -11) result = channel.bRatio();
    });
    if (result == 0.) print("no e+ e- channel for", particle.id());
    return result;
}

void MesonResonance::calcWidth(bool test) {
    if (debug) print("calcWidth", idRes, id1Abs, id2Abs, id3Abs, "with mass", particleDataPtr->m0(id1Abs));
    if (debug) print("multiplicity", mult);
    if (!has_neutrino(particlePtr->channel(iChannel))) { // SM channel, keeps its meMode
        widNow = GammaRes * particlePtr->channel(iChannel).bRatio();
        return;
    }
    auto id_lep = mult == 2 ? id2Abs : id3Abs;
    preFac = neutrino_coupling(id1Abs, id_lep + 1) * sqr(coupSMPtr->GF()) * Pythia8::pow3(mHat) / 8. / M_PI;
    if (debug) print("preFac", preFac, neutrino_coupling(id1Abs, id_lep + 1), sqr(coupSMPtr->GF()), Pythia8::pow3(mHat));
    widNow = 0.;
    if (preFac <= 0.) {
        if (debug) print("zero prefactor", idRes, id1Abs, id2Abs, id3Abs);
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
    switch (mult) {
        case 2 :
            if (is_heavy_neutral_lepton(id1Abs) && id2Abs > 10 && id2Abs < 17) {
                preFac *= (mr1 + mr2 - sqr(mr2 - mr1)) * std::sqrt(lambda(1., mr1, mr2));
                if (debug) print("preFac", preFac);
                if (idRes == 443 || idRes == 553) {
                    preFac *= 27. * mHat / 8. / M_PI / sqr(coupSMPtr->alphaEM(sqr(mHat)));
                    auto term = 4. / 3. * coupSMPtr->sin2thetaW();
                    preFac *= idRes == 443 ? sqr(1. - 2 * term) / 4. : sqr(1. - term);
                    widNow = preFac * electron_branching_ratio(*particlePtr) * particlePtr->mWidth();
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

NeutrinoResonance::NeutrinoResonance(Pythia8::Pythia& pythia, std::function<double (int id_heavy, int id_light)> const& coupling, double mass, int id, bool majorana_) :
    neutrino_coupling(coupling), majorana(majorana_) {
    if (debug) print("NeutrinoResonance", id);
    initBasic(id);
    particlePtr = pythia.particleData.particleDataEntryPtr(idRes);
    particlePtr->setM0(mass);
    particlePtr->setMMin(mass / 2);
    particlePtr->setMMax(mass * 2);
    if (debug) print("Particle", idRes, particlePtr->name(), "mass", particlePtr->m0(), "tau", particlePtr->tau0(), particlePtr->mWidth());
}

bool NeutrinoResonance::init(Pythia8::Info* infoPtrIn) {
    if (particlePtr->m0() >= MASSMIN) return ResonanceWidths::init(infoPtrIn);
    // the pointers that ResonanceWidths::init() would set before it rejects the mass
    infoPtr = infoPtrIn;
    settingsPtr = infoPtr->settingsPtr;
    particleDataPtr = infoPtr->particleDataPtr;
    loggerPtr = infoPtr->loggerPtr;
    coupSMPtr = infoPtr->coupSMPtr;
    coupSUSYPtr = infoPtr->coupSUSYPtr;
    calculate_widths();
    return true;
}

bool NeutrinoResonance::initBSM() {
    if (debug) print("initBSM");
    particlePtr->setIsResonance(false);
    particlePtr->setMayDecay(true);
    return true;
}

// HNL decays (1805.08567, section 3). Leptons as there. Hadrons (rewritten 24 Nov 2026):
// - one meson at all masses, with the single-meson widths of section 3.2: l- P+ (pi, K, D, D_s, B, B_c),
//   l- V+ (rho, K*, D*, D_s*, a_1), nu P0 (pi0, eta, eta', eta_c), nu V0 (rho0, omega, phi, J/psi, a_1^0);
// - for each quark current, the remainder of the inclusive width (quarks times 1 + Delta_QCD) over the one-meson
//   widths of that current, if positive and if the hadronic system can be heavier than the lightest states the
//   mesons above do not cover. The q qbar' pair is colour connected and fragments in Pythia's hadron level: meMode 94
//   (V-A, products in the order dbar u l-) for the charged current, 93 for the neutral one.
// At m_tau this reproduces tau decays: pi, rho and a_1 give 89% of the non-strange hadronic width, the remainder the rest.

namespace {

const double alpha_s_mZ = 0.118;

bool spin_one(int id) {
    return id % 10 == 3;
}

// the up- and down-type quarks of the charged current that makes a meson
std::pair<int, int> current_of(int meson) {
    switch (meson) {
        case 211 : case 213 : case 20213 : return {2, 1};
        case 321 : case 323 : return {2, 3};
        case 411 : case 413 : return {4, 1};
        case 431 : case 433 : return {4, 3};
        case 521 : return {2, 5};
        case 541 : return {4, 5};
        default : print("no charged current for meson", meson);
    }
    return {0, 0};
}

std::vector<int> charged_current_mesons(int up, int down) {
    for (auto meson : {211, 321, 411, 431, 521, 541}) if (current_of(meson) == std::make_pair(up, down)) {
        switch (meson) {
            case 211 : return {211, 213, 20213};
            case 321 : return {321, 323};
            case 411 : return {411, 413};
            case 431 : return {431, 433};
            default : return {meson};
        }
    }
    return {};
}

// neutral-current mesons of a quark flavour; u ubar and d dbar share the light ones (listed under d)
std::vector<int> neutral_current_mesons(int quark) {
    switch (quark) {
        case 1 : return {111, 221, 331, 113, 223, 20113};
        case 3 : return {333};
        case 4 : return {441, 443};
        default : return {};
    }
}

// Lightest hadronic mass [GeV] of the states that the mesons of a current do not cover. Light currents: tau spectral
// functions are saturated by pi, pi pi and 3 pi below about 1 GeV; s sbar by the phi below about 1.5 GeV; heavy
// currents: D pi, D K, B pi, B_c pi and D Dbar, B Bbar thresholds.
double remainder_threshold(int up, int down) {
    if (up == 4) return down == 1 ? 2.01 : down == 3 ? 2.36 : 6.42;
    return down == 5 ? 5.42 : 1.;
}

double remainder_threshold(int quark) {
    switch (quark) {
        case 1 : case 2 : return 1.;
        case 3 : return 1.5;
        case 4 : return 3.74;
        default : return 10.56;
    }
}

}

double L(double x2) {
    if (x2 < 0.000026112) return -21.; // The log becomes complex for electrons; this approximate value is irrelevant as it will be multiplied by a very small number.
    auto sqrt = std::sqrt(1. - 4. * x2);
    return std::log((1. - 3. * x2 - (1. - x2) * sqrt) / x2 / (1. + sqrt));
}

double NeutrinoResonance::coupling(int neutrino) {
    return neutrino_coupling(idRes, neutrino);
}

double NeutrinoResonance::qcd_factor() { // Delta_QCD of 1805.08567 with alpha_s(M_N) at three loops, frozen below m_tau where it is tested
    auto alpha = alpha_strong.alphaS(sqr(std::max(mHat, particleDataPtr->m0(15)))) / M_PI;
    return 1. + alpha + 5.2 * sqr(alpha) + 26.4 * cube(alpha);
}

double NeutrinoResonance::width_charged_meson(int lepton, int meson) {
    auto m_l = particleDataPtr->m0(lepton);
    auto m_h = particleDataPtr->m0(meson);
    if (mHat <= m_l + m_h) return 0.;
    auto x_l = sqr(m_l / mHat);
    auto x_h = sqr(m_h / mHat);
    auto [up, down] = current_of(meson);
    auto width = sqr(coupSMPtr->GF() * decay_constant(meson)) * CKM2(up, down) * Pythia8::pow3(mHat) / 16. / M_PI
        * std::sqrt(lambda(1., x_h, x_l)) * coupling(lepton + 1);
    if (spin_one(meson)) return width / sqr(m_h) * (sqr(1. - x_l) + x_h * (1. + x_l) - 2. * sqr(x_h));
    return width * (sqr(1. - x_l) - x_h * (1. + x_l));
}

double NeutrinoResonance::width_neutral_meson(int neutrino, int meson) {
    auto m_h = particleDataPtr->m0(meson);
    if (mHat <= m_h) return 0.;
    auto x_h = sqr(m_h / mHat);
    auto width = sqr(coupSMPtr->GF() * decay_constant(meson)) * Pythia8::pow3(mHat) / 32. / M_PI * sqr(1. - x_h)
        * coupling(neutrino);
    if (spin_one(meson)) return width * sqr(correction_factor(meson) / m_h) * (1. + 2. * x_h);
    return width;
}

// N -> l_alpha- U Dbar: quarks with N_c |V|^2 and the QCD factor, or leptons U = nu_beta, D = l_beta
double NeutrinoResonance::width_charged_fermions(int lepton, int up, int down) {
    if (mHat <= get_mass(lepton) + get_mass(up) + get_mass(down)) return 0.;
    auto width = sqr(coupSMPtr->GF()) * Pythia8::pow5(mHat) / 192. / Pythia8::pow3(M_PI) * coupling(lepton + 1)
        * three_body_width.get_width(idRes, lepton, up, down);
    return is_quark(up) ? width * 3. * CKM2(up, down) * qcd_factor() : width;
}

// N -> nu_alpha f fbar through the neutral current, and the charged current for f = l_alpha
double NeutrinoResonance::width_neutral_fermions(int neutrino, int fermion) {
    auto x = sqr(get_mass(fermion) / mHat);
    if (4. * x >= 1.) return 0.;
    auto root = std::sqrt(1. - 4. * x);
    auto term_1 = (1. - 14. * x - 2. * sqr(x) - 12. * cube(x)) * root + 12. * sqr(x) * (sqr(x) - 1.) * L(x);
    auto term_2 = x * (2. + 10. * x - 12 * sqr(x)) * root + 6 * sqr(x) * (1. - 2. * x + 2 * sqr(x)) * L(x);
    auto width = sqr(coupSMPtr->GF()) * Pythia8::pow5(mHat) / 192. / Pythia8::pow3(M_PI) * coupling(neutrino) * NZ(fermion)
        * (Cf1(neutrino, fermion) * term_1 + 4. * Cf2(neutrino, fermion) * term_2);
    return is_quark(fermion) ? width * qcd_factor() : width;
}

// N -> nu_alpha nu_beta nubar_beta, one of the three
double NeutrinoResonance::width_invisible(int neutrino, int neutrino_2) {
    return sqr(coupSMPtr->GF()) * Pythia8::pow5(mHat) / 192. / Pythia8::pow3(M_PI) * coupling(neutrino) / 4.
        * (neutrino == neutrino_2 ? 2. : 1.);
}

double NeutrinoResonance::remainder_charged(int lepton, int up, int down) {
    if (mHat - particleDataPtr->m0(lepton) < remainder_threshold(up, down)) return 0.;
    auto rest = width_charged_fermions(lepton, up, down);
    for (auto meson : charged_current_mesons(up, down)) rest -= width_charged_meson(lepton, meson);
    return std::max(rest, 0.);
}

double NeutrinoResonance::remainder_neutral(int neutrino, int quark) {
    if (mHat < remainder_threshold(quark)) return 0.;
    auto family = quark == 2 ? 1 : quark;
    auto inclusive = width_neutral_fermions(neutrino, family) + (family == 1 ? width_neutral_fermions(neutrino, 2) : 0.);
    auto rest = inclusive;
    for (auto meson : neutral_current_mesons(family)) rest -= width_neutral_meson(neutrino, meson);
    return rest > 0. ? rest * width_neutral_fermions(neutrino, quark) / inclusive : 0.;
}

bool NeutrinoResonance::allowCalc() {
    if (debug) print("allowCalc");
    mHat = particlePtr->m0(); // init() sets it only after allowCalc()
    particlePtr->clearChannels();
    for (auto lepton : charge_leptons()) {
        for (auto lepton_2 : charge_leptons()) if (lepton_2 != lepton) particlePtr->addChannel(1, 0., 94, -lepton_2, lepton_2 + 1, lepton);
        for (auto up : {2, 4}) for (auto down : {1, 3, 5}) {
            for (auto meson : charged_current_mesons(up, down)) particlePtr->addChannel(1, 0., 0, lepton, meson);
            particlePtr->addChannel(1, 0., 94, -down, up, lepton);
        }
    }
    for (auto neutrino : neutral_leptons()) {
        for (auto neutrino_2 : neutral_leptons()) particlePtr->addChannel(1, 0., 0, neutrino, neutrino_2, -neutrino_2);
        for (auto lepton : charge_leptons()) particlePtr->addChannel(1, 0., 0, neutrino, lepton, -lepton);
        for (auto quark : {1, 3, 4}) for (auto meson : neutral_current_mesons(quark)) particlePtr->addChannel(1, 0., 0, neutrino, meson);
        for (auto quark : {1, 2, 3, 4, 5}) particlePtr->addChannel(1, 0., 93, neutrino, quark, -quark);
    }
    if (majorana) add_charge_conjugates();
    return true;
}

void NeutrinoResonance::add_charge_conjugates() {
    auto size = particlePtr->sizeChannels();
    for (int pos = 0; pos < size; ++pos) {
        auto const& channel = particlePtr->channel(pos);
        std::array<int, 3> products{};
        for (int i = 0; i < channel.multiplicity(); ++i) {
            auto product = channel.product(i);
            products.at(i) = particleDataPtr->hasAnti(product) ? -product : product;
        }
        particlePtr->addChannel(1, 0., channel.meMode(), products[0], products[1], products[2]);
    }
}

void NeutrinoResonance::initConstants() {
    if (debug) print("initConstants");
    three_body_width.set_pointers(particleDataPtr, settingsPtr);
    alpha_strong.init(alpha_s_mZ, 3, 6, false);
}

double NeutrinoResonance::CKM2(int up, int down) {
    return coupSMPtr->V2CKMgen(up_idx(up), down_idx(down));
}

double NeutrinoResonance::Cf1(int neutrino, int fermion) {
    double first = is_quark(fermion) ? (is_up_type(fermion) ? 2. : 1.) / 3. : (fermion + 1 == neutrino ? -1 : 1);
    double second = is_quark(fermion) ? (is_up_type(fermion) ? 4. : 1.) / 9. : 1.;
    return 1. / 4. * (1. - 4. * first * coupSMPtr->sin2thetaW() + 8. * second * sqr(coupSMPtr->sin2thetaW()));
}

double NeutrinoResonance::Cf2(int neutrino, int fermion) {
    double first = is_quark(fermion) ? (is_up_type(fermion) ? 1. : 1. / 2.) / 3. : 1. / 2.;
    double second = is_quark(fermion) ? (is_up_type(fermion) ? 4. : 2.) / 3. : 2.;
    double third = is_quark(fermion) ? 1. : (fermion + 1 == neutrino ? -1. : 1.);
    return first * coupSMPtr->sin2thetaW() * (second * coupSMPtr->sin2thetaW() - third);
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
    auto const& channel = particlePtr->channel(iChannel);
    std::array<int, 3> ids{};
    for (int i = 0; i < channel.multiplicity(); ++i) ids.at(i) = std::abs(channel.product(i));
    if (channel.multiplicity() == 2) widNow = is_light_neutrino(ids[0]) ? width_neutral_meson(ids[0], ids[1]) : width_charged_meson(ids[0], ids[1]);
    else if (channel.meMode() == 94) widNow = is_quark(ids[1]) ? remainder_charged(ids[2], ids[1], ids[0]) : width_charged_fermions(ids[2], ids[1], ids[0]);
    else if (is_quark(ids[1])) widNow = remainder_neutral(ids[0], ids[1]);
    else if (is_light_neutrino(ids[1])) widNow = width_invisible(ids[0], ids[1]);
    else widNow = width_neutral_fermions(ids[0], ids[1]);
    if (debug) print("Calculated", widNow, "for", idRes, "to", ids[0], ids[1], ids[2]);
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
