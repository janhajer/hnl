#include "pythia.hh"
#include "ResonanceWidths.hh"
#include "branching_ratios.hh"

namespace hnl {

namespace {

constexpr bool debug = false;

}

void save_data(Result& result, hnl::Loop const& loop, double mass, int source) {
    std::ofstream output_file(std::to_string(source) + ".dat");
    output_file << 0 << '\t' << 1 << '\t' << 2 << '\t' << 3 << '\t' << 4;
    for (auto step = 0; step <= loop.steps; ++step) output_file << std::scientific << '\t' << loop.mass(mass, step);
    output_file << '\n';
    for (auto& row : result[source]) {
        output_file << std::scientific << std::get<0>(row.first) << '\t' << std::get<1>(row.first) << '\t' << std::get<2>(row.first) << '\t' << std::get<3>(row.first) << '\t' << std::get<4>(row.first) << '\t';
        for (auto step = 0; step <= loop.steps; ++step) output_file << std::scientific << row.second[step] << '\t';
        output_file << '\n';
    }
}

Result branching_ratio(Loop const& loop, double& mass_max, int source, int step) {
    print(step, "of", loop.steps);
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_branching_ratios(pythia);
    if (!is_heavy_neutral_lepton(source)) mass_max = pythia.particleData.m0(source);
    double coupling = 1;
    is_heavy_neutral_lepton(source) ? pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling(coupling), loop.mass(mass_max, step), source)) : pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling(coupling), source));
    pythia.particleData.particleDataEntryPtr(source)->rescaleBR();
    pythia.init();
    Result result;
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

void write_branching_ratios(int source) {
    Result result;
    Loop loop(.1, 5);
    double mass_max = 5;
    for (auto step = 0; step <= loop.steps; ++step) result += branching_ratio(loop, mass_max, source, step);
    save_data(result, loop, mass_max, source);
}

void write_branching_ratios() {

//     std::vector<int> sources{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541, 443, 553};
//     std::vector<int> sources{431, 411, 421};
//     std::vector<int> sources{511, 521, 531, 541};
//     std::vector<int> sources{511, 521, 531};
//     std::vector<int> sources{443, 553};
    std::vector<int> sources{heavy_neutrino};
    for (auto source : sources) write_branching_ratios(source);
}

}

