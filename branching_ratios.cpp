#include <map>
#include <array>

#include "pythia.hh"
#include "ResonanceWidths.hh"
#include "branching_ratios.hh"
#include "decay_table.hh"

namespace hnl {

namespace {

constexpr bool debug = false;

}

void save_data(BranchingRatios& result, hnl::Loop const& loop, double mass, int source) {
    std::ofstream output_file(std::to_string(source) + ".dat");
    output_file << 0 << '\t' << 1 << '\t' << 2 << '\t' << 3 << '\t' << 4;
    for (auto step = 0; step <= loop.steps; ++step) output_file << std::scientific << '\t' << loop.mass(mass, step);
    output_file << '\n';
    for (auto& row : result[source]) {
        output_file << std::scientific << row.first[0] << '\t' << row.first[1] << '\t' << row.first[2] << '\t' << row.first[3] << '\t' << row.first[4] << '\t';
        for (auto step = 0; step <= loop.steps; ++step) output_file << std::scientific << row.second[step] << '\t';
        output_file << '\n';
    }
}

void save_data(std::map<double,double>& result, hnl::Loop const& loop, double mass, int source) {
    std::ofstream output_file(std::to_string(source) + "lifetime.dat");
    output_file << "mass [GeV] " << '\t' << "Width [GeV]"  << '\n';
    for (auto& row : result) output_file << std::scientific << row.first << '\t' << row.second << '\n';
}

BranchingRatios branching_ratio(Loop const& loop, double& mass_max, int source, int step) {
    print("id", source, "step", step, "of", loop.steps);
//     print(meson_decay_table(neutrino_coupling(1), loop.mass(mass_max, step), heavy_neutrino, source));
//     return {};
//     return hnl_branching_ratios(neutrino_coupling(1), loop.mass(mass_max, step), source, step);
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_branching_ratios(pythia);
    if (!is_heavy_neutral_lepton(source)) mass_max = pythia.particleData.m0(source);
    double coupling = 1;


    auto mass = loop.mass(mass_max, step);
    set_pythia_stable(pythia,heavy_neutrino, mass);
    is_heavy_neutral_lepton(source) ? pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling(coupling), loop.mass(mass_max, step), source)) : pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling(coupling), source));

//     if(is_heavy_neutral_lepton(source)){
//         for (auto const& line : hnl_decay_table(neutrino_coupling(1), loop.mass(mass_max, step), source)) {
//         if(debug) print(line);
//         pythia.readString(line);
//     }
//     } else {
//             pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling(coupling), source));
//     }


    pythia.init();
    BranchingRatios result;
    auto const& particle = *pythia.particleData.particleDataEntryPtr(source);
    for_each(particle, [&result, source, step](Pythia8::DecayChannel const& channel){
        if (channel.bRatio() > 0. && (has_neutrino(channel) || is_heavy_neutral_lepton(source))) result[source][ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}][step] = channel.bRatio();
    });
//     for (auto pos = 0; pos < particle.sizeChannels(); ++pos) {
//         auto channel = particle.channel(pos);
//         auto ratio = channel.bRatio();
//         if (ratio > 0. && (has_neutrino(channel) || is_heavy_neutral_lepton(source)))
//             result[source][ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}][step] = ratio;
//     }
    return result;
}

std::map<double,double> lifetime(Loop const& loop, double& mass_max, int source, int step) {
    print("id", source, "step", step, "of", loop.steps);
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_branching_ratios(pythia);
    if (!is_heavy_neutral_lepton(source)) mass_max = pythia.particleData.m0(source);
    double coupling = 1;
    is_heavy_neutral_lepton(source) ? pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling(coupling), loop.mass(mass_max, step), source)) : pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling(coupling), source));
    pythia.init();
    std::map<double,double> result;
    auto const& particle = *pythia.particleData.particleDataEntryPtr(source);
    result[particle.m0()] = particle.mWidth();
    return result;
}

void write_branching_ratios(int source) {
    BranchingRatios result;
    Loop loop(.1, 20);
    double mass_max = 6.2;
    for (auto step = 0; step <= loop.steps; ++step) result += branching_ratio(loop, mass_max, source, step);
    save_data(result, loop, mass_max, source);
}

void write_branching_ratios() {
    std::vector<int> sources{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541, 443, 553};
//     std::vector<int> sources{211, 130, 310, 321};
//     std::vector<int> sources{431, 411, 421};
//     std::vector<int> sources{511, 521, 531, 541};
//     std::vector<int> sources{511, 521, 531};
//     std::vector<int> sources{443, 553};
//     std::vector<int> sources{heavy_neutrino};
    for (auto source : sources) write_branching_ratios(source);
}

void write_lifetime() {
    std::map<double,double> result;
    Loop loop(.1, 100);
    double mass_max = 6;
    for (auto step = 0; step <= loop.steps; ++step) result += lifetime(loop, mass_max, heavy_neutrino, step);
    save_data(result, loop, mass_max, heavy_neutrino);
}

}
