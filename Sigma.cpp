// Simple illustration how to provide (a) your own resonance-width class, and (b) your own cross-section class, with instances handed in to Pythia.
// The hypothetical scenario is that top would have been so long-lived that a toponium resonance Theta could form.
// Then production could proceed via q qbar -> gamma*/Z* -> Theta, with decay either to a fermion pair or (dominantly) to three gluons.
// The implementation is not physically correct in any number of ways, but should exemplify the strategy needed for realistic cases.

#include "Sigma.hh"

namespace hnl {

namespace {

const int debug = 0;

template <typename Value>
int sgn(Value value) {
    return (Value(0) < value) - (value < Value(0));
}

template<typename Object>
auto sqr(Object const& object) noexcept {
    return object * object;
}

template<typename Object>
auto cube(Object const& object) noexcept {
    return object * object * object;
}

void print() noexcept {
    std::cout << std::endl;
}

template<typename Object, typename ... Arguments>
void print(Object const& object, Arguments ... arguments) noexcept {
    std::cout << std::boolalpha << std::scientific << object << ' ';
    print(arguments ...);
}

}

Sigma::Sigma(int id_, double coupling_) {
    if (debug > 0) print("Sigma", id_);
    id = id_;
    coupling = coupling_;
}

int Sigma::code() const {
    if (debug > 1) print("code");
    return 99000;
}

std::string Sigma::inFlux() const {
    if (debug > 1) print("inFlux");
    return "qqbar";
}

void Sigma::initProc() {
    if (debug > 0) print("initProc", particleDataPtr->m0(id));
    mass2 = sqr(particleDataPtr->m0(id));
}

// int Sigma::resonanceA() const {
//     if (debug > 1) print("resonance A");
//     return 24;
// }

// int Sigma::resonanceB() const {
//     if (debug > 1) print("resonance B");
//     return 23;
// }

void Sigma::sigmaKin() { // should use Mandelstam and squares squared transvrse moemtum and nominal BreitWigner factors tH, uH, tH2, uH2, pT2, runBW3, runBW4;
    if (debug > 1) print("sigmaKin", sH, tH, uH, tH2, uH2, pT2);
    double ratio = mass2 / sH;
    if (ratio > 1) {
        if(debug > 1) print("too large mass");
        sigma = 0;
        return;
    } else {
        double kin = 1. - (3 * ratio - cube(ratio)) / 2;
        sigma = sH * sqr(couplingsPtr->GF()) * coupling / 6 / 3 / 4 / sqr(M_PI) * kin;
        if (debug > 1) print(sigma);
    }
}

double Sigma::sigmaHat() { // should use m3, s3, m4, s4;
    if (debug > 1) print("sigmaHat", id1, id2);
//     auto charge = particleDataPtr->chargeType(id1) + particleDataPtr->chargeType(id2);
//     auto lep = charge == 0 ? 12 : - sgn(charge) * 11;
//     auto openFracPair = particleDataPtr->resOpenFrac(id, lep);
//     print(openFracPair);
    return sigma * (id1 == -id2 ? 1. : couplingsPtr->V2CKMid(id1, id2));
}

void Sigma::setIdColAcol() {
    if (debug > 1) print("setIdColAcol", id1, id2);
    auto charge = particleDataPtr->chargeType(id1) + particleDataPtr->chargeType(id2);
    auto lep = charge == 0 ? 12 : - sgn(charge) * 11;
    setId(id1, id2, id, lep);
    setColAcol(1, 0, 0, 1, 0, 0, 0, 0); // Colour flow topologies.
    if (id1 < 0) swapColAcol(); // Swap when antiquarks.
}

double Sigma::weightDecay(Pythia8::Event& event, int iResBeg, int iResEnd) { // Evaluate weight for Theta -> g g g.
    return 1.;
}

std::string Sigma::name() const { // Info on the subprocess.
    if (debug > 1) print("name");
    return "qbar q' -> HNL lep";
}


bool Sigma::isSChannel() const {
    return true;
}

}

