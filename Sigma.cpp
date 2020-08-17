#include "Sigma.hh"
#include "math.hh"
#include "io.hh"

namespace hnl {

namespace {

const int debug = 1;

}

Sigma::Sigma(std::function<double (int id_heavy, int id_light)> const& neutrino_coupling_, int heavy_, int light_) :
    heavy(heavy_),
    light(light_),
    neutrino_coupling(neutrino_coupling_) {
    if (debug > 0) print("Sigma", heavy, light);
}

int Sigma::code() const {
    if (debug > 1) print("code");
    return 90000 + (heavy - 9900000) * 100 + light;
}

std::string Sigma::inFlux() const {
    if (debug > 0) print("inFlux");
    return "qqbar";
}

void Sigma::initProc() {
    if (debug > 0) print("initProc", particleDataPtr->m0(heavy));
    mass2 = sqr(particleDataPtr->m0(heavy));
}

void Sigma::sigmaKin() { // should use Mandelstam and squares squared transvrse moemtum and nominal BreitWigner factors tH, uH, tH2, uH2, pT2, runBW3, runBW4;
    if (debug > 1) print("sigmaKin", sH, tH, uH, tH2, uH2, pT2);
    double ratio = mass2 / sH;
    int Nc = 3;
    if (ratio > 1) {
        if (debug > 0) print("too large mass");
        M2 = 0;
    } else {
        M2 = 8. * Nc * sqr(couplingsPtr->GF()) * neutrino_coupling(heavy, light) * (uH - mass2);
        if (debug > 1) print(M2);
    }
}

double Sigma::sigmaHat() { // should use m3, s3, m4, s4;
    if (debug > 1) print("sigmaHat", id1, id2);
    auto charge = particleDataPtr->chargeType(id1) + particleDataPtr->chargeType(id2);
    auto lep = charge == 0 ? light : - sgn(charge) * (light - 1);
    auto resOpenFrac = particleDataPtr->resOpenFrac(heavy, lep);
    return M2 * (uH - sqr(particleDataPtr->m0(lep))) * (id1 == -id2 ? 1. : couplingsPtr->V2CKMid(id1, id2)) * resOpenFrac;
}

void Sigma::setIdColAcol() {
    if (debug > 1) print("setIdColAcol", id1, id2);
    auto charge = particleDataPtr->chargeType(id1) + particleDataPtr->chargeType(id2);
    auto lep = charge == 0 ? light : - sgn(charge) * (light - 1);
    setId(id1, id2, lep, heavy);
    setColAcol(1, 0, 0, 1, 0, 0, 0, 0);
    if (id1 < 0) swapColAcol();
}

double Sigma::weightDecay(Pythia8::Event&, int iResBeg, int iResEnd) {
//     if (debug > 0)
        print("weightDecay", iResBeg, iResEnd);
    return 1;
}

std::string switch_heavy(int heavy) {
    switch (heavy) {
        case 9900012 :
            return "N1";
        case 9900014 :
            return "N2";
        case 9900016 :
            return "N3";
        default :
            print("Sigma", "name", "do not end up here");
    }
    return "";
}

std::string switch_light(int light) {
    switch (light) {
        case 12 :
            return "e";
        case 14 :
            return "mu";
        case 16 :
            return "tau";
        default :
            print("Sigma", "name", "do not end up here");
    }
    return "";
}

std::string Sigma::name() const {
    if (debug > 1) print("name");
    return "qbar q' -> " + switch_light(light) + " " + switch_heavy(heavy);
}

int Sigma::id3Mass() const {
    if (debug > 1) print("id3Mass", id1, id2);
    if (id1  == 0 || id2 == 0) { // BUG
        if (debug > 0) print("wrong mass");
        return light - 1;
    }
    auto charge = particleDataPtr->chargeType(id1) + particleDataPtr->chargeType(id2);
    return charge == 0 ? light : - sgn(charge) * (light - 1);
}

int Sigma::id4Mass() const {
    if (debug > 0) print("id4Mass", id1, id2);
    return heavy;
}

}

