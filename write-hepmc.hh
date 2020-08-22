#pragma once

#include "ResonanceWidths.hh"
#include "Sigma.hh"
#include "hepmc.hh"
#include "string.hh"
#include "pythia.hh"

namespace hnl {

namespace {

const bool debug = false;

}

std::pair<int, double> get_max_width(Pythia8::Pythia& pythia, std::vector<int> const& mesons) {
    double max_width = 0;
    int id;
    for (auto meson : mesons) {
        double partial_width = 0.;
        auto& particle = *pythia.particleData.particleDataEntryPtr(meson);
       for_each(particle, [&partial_width, meson](Pythia8::DecayChannel const& channel){
            if (has_neutrino(channel) && channel.bRatio() > 0.) {
                partial_width += channel.bRatio();
                if (debug) print(meson, channel.product(0), channel.product(1), channel.product(2));
            }
        });
//         for (auto number : irange(particle.sizeChannels())) {
//             auto& channel = particle.channel(number);
//             if (has_neutrino(channel) && channel.bRatio() > 0.) {
//                 partial_width += channel.bRatio();
//                 if (debug) print(meson, channel.product(0), channel.product(1), channel.product(2));
//             }
//         }
        if (max_width < partial_width) {
            max_width = partial_width;
            id = meson;
        }
    }
    return {id, max_width};
}

auto get_optimal_coupling(double mass, std::vector<int> const& mesons) {
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_write_hepmc(pythia,heavy_neutrino, mass);
    set_pythia_init_quiet(pythia);
    for (auto meson : mesons) pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling(1), meson));
    pythia.init();

    auto [id, max_width] = get_max_width(pythia, mesons);

    print("The maximal BR into HNLs with U^2 = 1 is", max_width, "from", id);
    if (mass < 0.13957) return 1. / max_width / 10000;
    if (mass < 0.49761) return 1. / max_width / 100;
    if (mass > 5) return 1. / max_width * 10000;
    if (mass > 4) return 1. / max_width * 1000;
    if (mass > 3) return 1. / max_width * 100;
    if (mass > 2.55) return 1. / max_width * 10;
    return 1. / max_width;
}

void for_each_if(Pythia8::Pythia& pythia, HepMC::IO_GenEvent& file, std::function<bool(Pythia8::Event const&)> const& function) {
    HepMC::Pythia8ToHepMC converter;
    int event_number = 0;
    while (event_number < pythia.mode("Main:numberOfEvents")) {
        if (!pythia.next()) continue;
        if (!function(pythia.event)) continue;
        ++event_number;
        HepMC::GenEvent event;
        converter.fill_next_event(pythia, &event);
        file.write_event(&event);
    }
}

void write_hepmc(Pythia8::Pythia& pythia, double mass, double coupling) {
    HepMC::IO_GenEvent file(std::to_string(mass) + ".hep", std::ios::out);
    file.write_comment("mass " + std::to_string(mass) + " GeV");

    for (auto heavy : heavy_neutral_leptons()) for (auto light : light_neutrinos()) file.write_comment("coupling " + std::to_string(heavy) + " " + std::to_string(light) + " " + std::to_string(neutrino_coupling(coupling)(heavy, light)));

    int total = 0;
    int successfull = 0;
    int too_many = 0;

    for_each_if(pythia, file, [&too_many, &total, &successfull](Pythia8::Event const & event) -> bool {
        ++total;
        int found_one = 0;
        bool success = false;
        for (auto line = 0; line < event.size(); ++line) {
            auto const& particle = event[line];
            if (!is_heavy_neutral_lepton(std::abs(particle.id()))) continue;
            if (!is_heavy_neutral_lepton(std::abs(event[particle.mother1()].id()))) ++found_one;
            if (found_one > 1) {
                print("More than one neutrino in event", total, "in line", line, "from mother", event[particle.mother1()].id());
                success = false;
                break;
            }
            double eta = std::abs(particle.eta());
            if (!(1.4 < eta && eta < 3.5)) continue;
            if (!is_heavy_neutral_lepton(event[particle.mother1()].id())) print("Success event", successfull + 1, "of", total, "in line", line, "from mother", event[particle.mother1()].id(), "with decay vertex", particle.vDec().pAbs(), "and eta", particle.eta(), "and phi", particle.phi());
            success = true;
        }
        if (found_one > 1) ++too_many;
        if (success) ++successfull;
        return success;
    });

    pythia.stat();

    print("Successfully saved", successfull, "events");
    print(successfull + too_many, "of", total, "events had at least one HNL. That is a fraction of", double(successfull + too_many) / total);
    print("Therefore of", pythia.info.sigmaGen(), "mb we save", pythia.info.sigmaGen() * (successfull + too_many) / total, "mb");

    file.write_comment("sigma " + to_string(pythia.info.sigmaGen() * successfull / total) + " mb");
}

void write_hepmc(double mass) {
    print("Generating events for HNLs with mass", mass, "GeV");

    std::vector<int> mesons{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541, 443, 553};
    auto coupling = get_optimal_coupling(mass, mesons);
    print("Use U^2 =", coupling);

    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_write_hepmc(pythia,heavy_neutrino, mass);
    for (auto meson : mesons) pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling(coupling), meson));

    pythia.init();

    auto [id, max_width] = get_max_width(pythia, mesons);
    print("The maximal BR into HNLs with U^2 =", coupling, "is", max_width, "from", id);

    write_hepmc(pythia, mass, coupling);
}

void write_hepmcs() {
    Loop loop(.1, 100);
    for (auto step = 0; step <= loop.steps; ++step) {
        auto mass = loop.mass(6., step);
        std::ofstream ofstream(std::to_string(mass) + ".txt");
        std::streambuf* streambuf = std::cout.rdbuf();
        std::cout.rdbuf(ofstream.rdbuf());
        write_hepmc(mass);
        std::cout.rdbuf(streambuf);
    }
}

void calculate_sigma(double mass) {
    if (debug) print("calculate_sigma");
    auto coupling = 1;
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_sigma(pythia);
    pythia.readString("Main:numberOfEvents = 100");
    pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling(coupling), mass, heavy_neutrino));
    auto sigma = std::make_unique<Sigma>(neutrino_coupling(coupling), heavy_neutrino, 12);
    pythia.setSigmaPtr(sigma.get());
    pythia.init();
    for (auto event_number = 0; event_number < pythia.mode("Main:numberOfEvents"); ++event_number) if (!pythia.next()) print("Error in event", event_number);
    pythia.stat();
    print("sigma", pythia.info.sigmaGen() / coupling);
}

void write_sigma_hepmc(double mass) {
    print("Generating events for HNLs with mass", mass, "GeV");

    auto coupling = 1;
    print("Use U^2 =", coupling);

    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_sigma(pythia);
    pythia.readString("Main:numberOfEvents = 100000");

    pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling(coupling), mass, heavy_neutrino));
    auto sigma = std::make_unique<Sigma>(neutrino_coupling(coupling), heavy_neutrino, 12);
    pythia.setSigmaPtr(sigma.get());

    pythia.init();

    write_hepmc(pythia, mass, coupling);

}

void write_sigma_hepmcs() {
    Loop loop(.1, 10);
    for (auto step = 0; step <= loop.steps; ++step) {
        auto mass = loop.mass(6.2, step);
        std::ofstream ofstream("sigma_" + std::to_string(mass) + ".txt");
        std::streambuf* streambuf = std::cout.rdbuf();
        std::cout.rdbuf(ofstream.rdbuf());
        write_sigma_hepmc(mass);
        std::cout.rdbuf(streambuf);
    }
}

auto get_optimal_coupling_2(double mass, std::vector<int> const& mesons) {
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_minimum_bias(pythia,heavy_neutrino, mass);
    set_pythia_init_quiet(pythia);
    for (auto meson : mesons) pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling(1), meson));
    pythia.init();

    auto [id, max_width] = get_max_width(pythia, mesons);

    print("The maximal BR into HNLs with U^2 = 1 is", max_width, "from", id);
    if (mass < 0.13957) return 1. / max_width / 10000;
    if (mass < 0.49761) return 1. / max_width / 100;
    if (mass > 5) return 1. / max_width * 10000;
    if (mass > 4) return 1. / max_width * 1000;
    if (mass > 3) return 1. / max_width * 100;
    if (mass > 2.3) return 1. / max_width * 10;
    return 1. / max_width;
}

void write_minimum_bias(double mass) {
    print("Generating events for HNLs with mass", mass, "GeV");

    std::vector<int> mesons{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541, 443, 553};
    auto coupling = get_optimal_coupling_2(mass, mesons);
    print("Use U^2 =", coupling);

    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_minimum_bias(pythia,heavy_neutrino, mass);
    for (auto meson : mesons) pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling(coupling), meson));

    pythia.init();

    auto [id, max_width] = get_max_width(pythia, mesons);
    print("The maximal BR into HNLs with U^2 =", coupling, "is", max_width, "from", id);

    write_hepmc(pythia, mass, coupling);
}

void write_minimum_biases() {
    Loop loop(.1, 100);
    for (auto step = 0; step <= loop.steps; ++step) {
        auto mass = loop.mass(6., step);
        std::ofstream ofstream(std::to_string(mass) + ".txt");
        std::streambuf* streambuf = std::cout.rdbuf();
        std::cout.rdbuf(ofstream.rdbuf());
        write_minimum_bias(mass);
        std::cout.rdbuf(streambuf);
    }
}

}

