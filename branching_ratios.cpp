#include <array>
#include <fstream>
#include <map>

#include "pythia.hh"
#include "ResonanceWidths.hh"
#include "branching_ratios.hh"

namespace hnl {

namespace {

// source -> decay products -> mass step -> branching ratio
using BranchingRatios = std::map<int, std::map<std::array<int, 5>, std::map<int, double>>>;

void save_data(BranchingRatios& result, Loop const& loop, double mass, int source) {
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

void save_data(std::map<double, double> const& result, int source) {
    std::ofstream output_file(std::to_string(source) + "lifetime.dat");
    output_file << "mass [GeV] " << '\t' << "Width [GeV]" << '\n';
    for (auto const& row : result) output_file << std::scientific << row.first << '\t' << row.second << '\n';
}

// HNL decays, or meson decays into HNLs of all three flavours with the same mass
void add_resonance(Pythia8::Pythia& pythia, int source, double mass) {
    double coupling = 1;
    if (is_heavy_neutral_lepton(source)) {
        pythia.addResonancePtr(std::make_shared<NeutrinoResonance>(pythia, neutrino_coupling(coupling), mass, source));
        return;
    }
    for (auto hnl : heavy_neutral_leptons()) set_pythia_stable(pythia, hnl, mass);
    pythia.addResonancePtr(std::make_shared<MesonResonance>(pythia, neutrino_coupling(coupling), source));
}

BranchingRatios branching_ratio(Loop const& loop, double& mass_max, int source, int step) {
    print("id", source, "step", step, "of", loop.steps);
    Pythia8::Pythia pythia(xml_dir(), false);
    set_pythia_branching_ratios(pythia);
    if (!is_heavy_neutral_lepton(source)) mass_max = pythia.particleData.m0(source);
    add_resonance(pythia, source, loop.mass(mass_max, step));
    pythia.init();
    BranchingRatios result;
    auto const& particle = *pythia.particleData.particleDataEntryPtr(source);
    for_each(particle, [&result, source, step](Pythia8::DecayChannel const& channel) {
        if (channel.bRatio() > 0. && (has_neutrino(channel) || is_heavy_neutral_lepton(source))) result[source][ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}][step] = channel.bRatio();
    });
    return result;
}

std::map<double, double> lifetime(Loop const& loop, double mass_max, int source, int step) {
    print("id", source, "step", step, "of", loop.steps);
    Pythia8::Pythia pythia(xml_dir(), false);
    set_pythia_branching_ratios(pythia);
    add_resonance(pythia, source, loop.mass(mass_max, step));
    pythia.init();
    auto const& particle = *pythia.particleData.particleDataEntryPtr(source);
    return {{particle.m0(), particle.mWidth()}};
}

}

void write_branching_ratios(int source) {
    BranchingRatios result;
    Loop loop(.1, 100);
    double mass_max = 6.2; // meson sources: their mass
    for (auto step = 0; step <= loop.steps; ++step) result += branching_ratio(loop, mass_max, source, step);
    save_data(result, loop, mass_max, source);
}

void write_branching_ratios() {
    // other sources: 211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541, 443, 553
    std::vector<int> sources{heavy_neutral_electron};
    for (auto source : sources) write_branching_ratios(source);
}

void write_lifetime() {
    std::map<double, double> result;
    Loop loop(.1, 100);
    for (auto step = 0; step <= loop.steps; ++step) result += lifetime(loop, 6., heavy_neutral_lepton, step);
    save_data(result, heavy_neutral_lepton);
}

}
