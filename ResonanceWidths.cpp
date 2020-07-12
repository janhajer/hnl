// #include "generic.hh"
#include "ResonanceWidths.hh"

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
    std::cout << std::boolalpha << std::scientific << object << ' ';
    print(arguments ...);
}

}

bool MesonResonance::initBSM()
{
    print("initBSM");
//     standard_model = dynamic_cast<Pythia8::CoupSM*>(couplingsPtr);
    return true;
}

bool MesonResonance::allowCalc()
{
    print("allowCalc");
    bool done = getChannels(idRes);
    std::stringstream idStream;
    idStream << "ID = " << idRes ;
    if (!done) infoPtr->errorMsg("Error in SusyResonanceWidths::allowcalc: " "unable to reset decay table.", idStream.str(), true);
    return done;
}

bool MesonResonance::getChannels(int resid)
{
    print("getChannels");
    id = std::abs(resid);
    Pythia8::ParticleDataEntry* particle_data_entry = particleDataPtr->particleDataEntryPtr(id);
//     particle_data_entry->clearChannels(); // Delete any decay channels read
    std::vector<int> neutrinos{14};
    for (auto const& neutrino : neutrinos) {
        if (particle_data_entry->chargeType() == 0) {
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino, 111);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino - 1, 211);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino, 311);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino - 1, 321);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino, 411);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino - 1, 421);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino, 113);
        } else {
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino - 1);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino - 1, 111);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino, 211);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino - 1, 311);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino, 321);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino - 1, 411);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino, 421);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino, 213);
            particle_data_entry->addChannel(1, 0., 0, 9900014, neutrino, 413);
        }
    }
    return true;
}

void MesonResonance::initConstants() // Initialize constants.
{
    print("initConstants");
    pseudo_scalar.setPointers(particleDataPtr, &standard_model, infoPtr); // Initialize functions for calculating 3-body widths
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
        return {0, 0};
    }
}

double MesonResonance::NeutU(int id_heavy, int id_light)
{
    auto diff = id_heavy - id_light - 9900000;
    print("diff", diff);
    double cpl = (diff == 0 || diff == 1) ? neutrino_coupling : 0.;
    print("coupl", cpl);
    return cpl;
}

double clepsch_gordan_2(int id)
{
    return id == 221 ? .5 : 1;
}

double MesonResonance::CKM2(int id)
{
    std::pair<int, int> pair = quark_pair(id);
//     return .1;
    return standard_model.V2CKMgen(pair.first, pair.second);
}

void MesonResonance::calcPreFac(bool) // Common coupling factors.
{
    print("calcPreFac");
}

double lambda(double a, double b, double c)
{
    return sqr(a) + sqr(b) + sqr(c) - 2 * a * b - 2 * a * c - 2 * b * c;
}

void MesonResonance::calcWidth(bool test) // Calculate width for currently considered channel.
{
//     if(test) return;
    print("calcWidth", id1Abs, id2Abs);
    double GF2 = sqr(standard_model.GF());
    if (ps == 0.) return; // Check that mass is above threshold.
    print("phasespace");
    if (mult == 2) { // Two-body decays
        if (id1Abs == 9900014 && id2Abs > 10 && id2Abs < 17) {
            preFac = GF2 * sqr(decay_constant(id)) * Pythia8::pow3(mHat) / 8 / M_PI * CKM2(id) * NeutU(id1, id2);
            widNow = preFac * (mr1 + mr2 - sqr(mr2 - mr1)) * std::sqrt(lambda(1, mr1, mr2));
        } else {
            print("Two-body not implemented");
        }
    } else {
        if (id1Abs == 9900014 && id3Abs > 10 && id3Abs < 17 && id2Abs > 20) {
            if(id2Abs == 113 || id2Abs == 213 || id2Abs == 313 || id2Abs == 323 || id2Abs == 413 || id2Abs == 423 || id2Abs == 433){
            preFac = GF2 * Pythia8::pow7(mHat) / 64 / cube(M_PI) / sqr(mf2) * clepsch_gordan_2(id2Abs) * CKM2(id) * NeutU(id1, id3);
            }else {
            preFac = GF2 * Pythia8::pow5(mHat) / 64 / cube(M_PI) * clepsch_gordan_2(id2Abs) * CKM2(id) * NeutU(id1, id3);
            }
            widNow = preFac * pseudo_scalar.getWidth(idRes, id2Abs, id3Abs);
        } else {
            print("Three-body for", id1Abs, id2Abs, id3Abs, "not implemented");
        }
    }
}

}
