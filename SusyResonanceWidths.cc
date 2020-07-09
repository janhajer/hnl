#include "SusyResonanceWidths.h"
#include "generic.hh"

namespace neutrino
{

bool SUSYResonanceWidths::initBSM()
{
    standard_model = (Pythia8::CoupSM*) couplingsPtr;
    return true;
}

bool SUSYResonanceWidths::allowCalc()
{
    bool done = getChannels(idRes);
    std::stringstream idStream;
    idStream << "ID = " << idRes ;
    if (!done) infoPtr->errorMsg("Error in SusyResonanceWidths::allowcalc: " "unable to reset decay table.", idStream.str(), true);
    return done;
}

bool ResonanceSlepton::getChannels(int resid)
{
    id = std::abs(resid);
    Pythia8::ParticleDataEntry* particle_data_entry = particleDataPtr->particleDataEntryPtr(id);
    particle_data_entry->clearChannels(); // Delete any decay channels read
    if (id % 2 == 1) {
        particle_data_entry->addChannel(1, 0., 0, 9900014, 11);
        particle_data_entry->addChannel(1, 0., 0, 9900014, 13);
        particle_data_entry->addChannel(1, 0., 0, 9900014, 15);
    } else {
        particle_data_entry->addChannel(1, 0., 0, 9900014, 12);
        particle_data_entry->addChannel(1, 0., 0, 9900014, 14);
        particle_data_entry->addChannel(1, 0., 0, 9900014, 16);
    }
    return true;
}

void ResonanceSlepton::initConstants() // Initialize constants.
{
    gf2 = sqr(standard_model->GF());
    stauWidths.setPointers(particleDataPtr, standard_model, infoPtr); // Initialize functions for calculating 3-body widths
}

double decay_constant(int id)
{
    switch (id) {
    case 211 : return 130.2;
    case 321 : return 155.6;
    case 411 : return 212.;
    case 413 : return 249.;
    case 521 : return 187.;
    case 541 : return 434.;
    case 111 : return 130.2;
    case 221 : return 81.7;
    case 331 : return -94.7;
    case 441 : return 237.;
    default : print("meson ", id, " not known");
    return 0.;
    }
}

std::pair<int, int> quark_pair(int id)
{
    switch (id) {
    case 211 : return {1, 1};
    case 321 : return {1, 2};
    case 411 : return {2, 1};
    case 413 : return {2, 2};
    case 521 : return {1, 3};
    case 541 : return {2, 3};
    case 111 : return {1, 1};
    case 221 : return {1, 2};
    case 331 : return {1, 2};
    case 441 : return {2, 2};
    default : print("meson ", id, " not known");
    return {0,0};
    }
}

double NeutU(int id){
    return id == 13 ? 1e-6 : 0.;
}

double clepsch_gordan(id){
    return id == 221 ? M_SQRT1_2 : 1;
}

double CKM2(int id)
{
    std::pair<int, int> pair = quark_pair(id);
    return standard_model->V2CKMgen(pair.first, pair.second);
}

void ResonanceSlepton::calcPreFac(bool) // Common coupling factors.
{
    alpEM = standard_model->alphaEM(mHat * mHat);
    preFac = gf2 * sqr(decay_constant(id)) * Pythia8::pow3(mHat) / 8 / M_PI * CKM2(id) * NeutU(id2);
}

double lambda(double a, double b, double c)
{
    return sqr(a) + sqr(b) + sqr(c) - 2 * a * b - 2 * a * c - 2 * b * c;
}

double Lambda(double xi)
{
    return std::sqrt(lambda(1, sqr(), xi))
}

void ResonanceSlepton::calcWidth(bool) // Calculate width for currently considered channel.
{
    if (ps == 0.) return; // Check that mass is above threshold.

    if (mult == 2) { // Two-body decays
        if (id1Abs <  && id2Abs < 17) {
            widNow = preFac * (mr2 + mr1 - sqr(mr2 - mr1)) * std::sqrt(lambda(1, mr2, mr1));
        }
    } else {
        preFac = gf2 * Pythia8::pow5(mHat) / 64 / cube(M_PI) * clepsch_gordan(id2Abs) * CKM2(id) * NeutU(id2);
        widNow = preFac * stauWidths.getWidth(idRes, id3Abs);
    }
}

}

