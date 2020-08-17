#pragma once

#include "hepmc.hh"
#include "geometry.hh"
#include "ResonanceWidths.hh"
#include "pythia.hh"
#include "read-file.hh"
#include "string.hh"

namespace hnl {

auto max(std::map<int, std::map<int, double>> const& couplings) {
    double max = 0.;
    for (auto const& inner : couplings) for (auto const& pair : inner.second) if (pair.second > max) max = pair.second;
    return max;
}

auto to_cgal(Pythia8::Particle const& particle) -> cgal::Point {
    return {particle.xProd() / 1000, particle.yProd() / 1000, particle.zProd() / 1000}; // convert from mm to m
}

void for_each_until(HepMC::GenEvent const& event, std::function<bool(HepMC::GenParticle const&)> const& function) {
    for (auto particle = event.particles_begin(); particle != event.particles_end(); ++particle) if (function(**particle)) return;
    print("no neutrino found");
}

Pythia8::Particle retrive_neutrino(HepMC::GenEvent const& event, double lifetime) {
    Pythia8::Particle pythia_particle;
    for_each_until(event, [&pythia_particle](HepMC::GenParticle const & hep_particle) {
        if (!is_heavy_neutral_lepton(hep_particle.pdg_id())) return false;
        auto& momentum = hep_particle.momentum();
        pythia_particle = Pythia8::Particle(hep_particle.pdg_id(), hep_particle.status(), 0, 0, 0, 0, 0, 0, to_pythia(momentum), momentum.m(), momentum.m(), 9.);
        pythia_particle.vProd(to_pythia(hep_particle.production_vertex()->position()));
        return true;
    });
    pythia_particle.tau(lifetime);
    return pythia_particle;
}

void for_each_until(HepMC::IO_GenEvent& events, std::function<bool(HepMC::GenEvent const&)> const& function) {
    auto* event = events.read_next_event();
    if (!event) print("Hepmc file is empty");
    while (event) {
        if (function(*event)) break;
        delete event;
        events >> event;
    }
}

double read_hepmc(boost::filesystem::path const& path, Meta const& meta, double coupling) {
    if (debug) print("read hep mc", path.string(), "with", coupling);
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_read_hepmc(pythia);

    auto couplings = [&meta, coupling](int heavy, int light) {
        return meta.couplings.at(heavy).at(light) > 0 ? coupling : 0.;
    };
    pythia.setResonancePtr(new NeutrinoResonance(pythia, couplings, meta.mass, heavy_neutrino));
    pythia.init();

    auto lifetime = pythia.particleData.tau0(heavy_neutrino);
    if (debug) print("trying to open", path.string());
    HepMC::IO_GenEvent events(path.string(), std::ios::in);
    if (debug) print("with result", events.error_message());
    int total = 0;
    int good = 0;
    int events_max = 1e6;
    auto analysis = mapp::analysis();
    for_each_until(events, [&](HepMC::GenEvent const & event) -> bool {
        ++total;
        pythia.event.reset();
        pythia.event.append(retrive_neutrino(event, lifetime));
        if (!pythia.next()) {
            print("Pythia encountered a problem");
            return false;
        }
        if (debug) pythia.event.list(true);
        for (auto line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            if (particle.chargeType() == 0) continue;
            if (!analysis.is_inside(to_cgal(particle))) continue;
            if (debug) print("Hooray!", particle.vProd().pAbs());
            ++good;
            break;
        }
        return total > events_max;
    });
    print("HNLs with m =", meta.mass, "GeV");
    auto result = meta.sigma;
    print("The production cross section is", result, "mb");
    auto rescaling = coupling / max(meta.couplings);
    print("Events were generated with U^2 =", max(meta.couplings), "Decays are porformed for U^2 =", coupling, "Hence the cross section is rescaled by", rescaling);
    result = meta.sigma * rescaling;
    print("The rescaled cross section is", result, "mb");
    auto fraction = double(good) / total;
    print(good, "of", total, "events accepted, that is a fraction", fraction);
    result *= fraction;
    print("In super MAPP", result, "mb");
    result /= 16.;
    print("In MAPP", result, "mb");
    return result;
}

using ScanResult = std::map<double, std::map<double, double>>;

void save_result(ScanResult const& result, std::string const& name = "result") {
    std::ofstream file;
    file.open(name + ".dat");
    bool first = true;
    for (auto const& line : result) {
        if (first) {
            file << "mass" << '\t';
            for (auto const& cell : line.second) file << cell.first << '\t';
            file << std::endl;
            first = false;
        }
        file << line.first << '\t';
        for (auto const& cell : line.second) file << cell.second << '\t';
        file << std::endl;
    }
}

double read_hepmc(boost::filesystem::path const& path, double coupling) {
    auto meta = meta_info(path);
    return meta ? read_hepmc(path, *meta, coupling) : 0.;
}

ScanResult scan_hepmc(boost::filesystem::path const& path) {
    ScanResult result;
    print("file", path);
    auto meta = meta_info(path);
    if (meta) for (auto coupling : log_range(1e-8, 1, 80)) result[meta->mass][coupling] = read_hepmc(path, *meta, coupling);
    return result;
}

void scan_hepmc(std::string const& path_name) {
    boost::filesystem::path path("./" + path_name);
    save_result(scan_hepmc(path), path.stem().string());
}

void scan_hepmcs(std::string const& path_name) {
    if (debug) print("read hep mcs in", path_name);
    ScanResult result;
    for (auto const& file : boost::make_iterator_range(boost::filesystem::directory_iterator(path_name), {})) if (file.path().extension().string() == ".hep") result += scan_hepmc(file.path());
    save_result(result);
}

}

