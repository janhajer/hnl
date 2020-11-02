#include "pythia.hh"

#include "decay_table.hh"
#include "ResonanceWidths.hh"
#include "string.hh"

namespace hnl {

namespace {

bool debug = false;

using Result = std::map<int, std::map<std::array<int, 5>, std::map<int, double>>>;

}

std::vector<std::string> decay_table(Pythia8::ParticleDataEntry const& particle, int meMode = -1) {
    if (debug) print("decay table");
    std::vector<std::string> result({std::to_string(particle.id()) + ":new = " + particle.name(particle.id()) + " " + particle.name(-particle.id()) + " " + std::to_string(particle.spinType()) + " " + std::to_string(particle.chargeType()) + " " + std::to_string(particle.colType()) + " " + std::to_string(particle.m0()) + " " + to_string(particle.mWidth()) + " " + to_string(particle.mMin()) + " " + std::to_string(particle.mMax()) + " " + to_string(particle.tau0())});
    for_each(particle, [&result, &particle, meMode](Pythia8::DecayChannel const & channel) {
        std::string string = std::to_string(particle.id()) + ":addChannel = " + std::to_string(channel.onMode()) + " " + std::to_string(channel.bRatio()) + " " + std::to_string(meMode > -1 ? meMode : channel.meMode());
        for_each(channel, [&string](int product) {
            string += " " + std::to_string(product);
        });
        result.emplace_back(string);
    });
    return result;
}

std::vector< std::string > meson_decay_table(std::function< double (int heavy, int light) > const& coupling, double mass, int id, int meson_id) {
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_stable(pythia, id, mass);
    set_pythia_passive(pythia);
    auto resonance = MesonResonance(pythia, coupling, meson_id);
    pythia.init();
    return decay_table(resonance.calculate_widths());
}

std::vector<std::string> hnl_decay_table(std::function<double (int heavy, int light) > const& coupling, double mass, int id) {
    if (debug) print("hnl decay table");
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_passive(pythia);
    pythia.init();
    auto resonance = NeutrinoResonance(pythia, coupling, mass, id);
    return decay_table(resonance.calculate_widths(), 101);
}

BranchingRatios hnl_branching_ratios(std::function<double (int heavy, int light) > const& coupling, double mass, int id, int step) {
    if (debug) print("hnl branching ratios", mass);
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_passive(pythia);
    if (!is_heavy_neutral_lepton(id)) set_pythia_stable(pythia, heavy_neutral_lepton, mass);
    pythia.init();
    BranchingRatios branching_ratios;
    print("Init mes res");
    auto res = MesonResonance(pythia, coupling, id);
    print("got res starting calculate_widths");
    auto& particle =
//     is_heavy_neutral_lepton(id) ? NeutrinoResonance(pythia, coupling, mass, id).calculate_widths() :
    res.calculate_widths();
    for_each(particle, [&branching_ratios, id, step](Pythia8::DecayChannel const & channel) {
//         if (channel.bRatio() > 0.
//             && (has_neutrino(channel) || is_heavy_neutral_lepton(id))
//         )
        branching_ratios[id][ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}][step] = channel.bRatio();
    });
    if (debug) print("BRs:");
    if (debug) print(branching_ratios);
    return branching_ratios;
}

}
