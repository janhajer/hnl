#pragma once

#include "read-file.hh"
#include "pythia-cgal.hh"
#include "mapp.hh"
#include "ResonanceWidths.hh"
#include "hepmc.hh"
#include "string.hh"

namespace hnl {

namespace hepmc {

namespace {

const bool debug = false;

}

auto sigma(std::vector<std::string> const& lines) {
    if (debug) print("find sigma");
    return find_if(lines, 1, [](auto const & strings)  {
        return /*strings.size() == 3 &&*/ strings.at(0) == "sigma" && strings.at(2) == "mb";
    });
}

auto mass(std::vector<std::string> const& lines) {
    if (debug) print("find mass");
    return find_if(lines, 1, [](auto const & strings)  {
        return strings.size() == 3 && strings.at(0) == "mass" && strings.at(2) == "GeV";
    });
}

auto coupling(std::vector<std::string> const& lines, int heavy, int light) {
    if (debug) print("find coupling");
    return find_if(lines, 3, [&](auto const & strings)  {
        return strings.size() == 4 && strings.at(0) == "coupling" && strings.at(1) == std::to_string(heavy) && strings.at(2) == std::to_string(light);
    });
}


boost::optional<Meta> meta_info(boost::filesystem::path const& path) {
    auto lines = import_head(path, 100) + import_tail(path, 100);
    Meta meta;
    meta.mass = to_double(mass(lines));
    if (meta.mass <= 0) return boost::none;
    meta.sigma = to_double(sigma(lines));
    if (meta.sigma <= 0) return boost::none;
    for (auto heavy : heavy_neutral_leptons()) for (auto light : light_neutrinos()) meta.couplings[heavy][light] = to_double(coupling(lines, heavy, light));
    if (meta.couplings.empty()) return boost::none;
    print("Meta info", meta.mass, meta.sigma);
    return meta;
}

void for_each_until(HepMC::GenEvent const& event, std::function<bool(HepMC::GenParticle const&)> const& function) {
    for (auto particle = event.particles_begin(); particle != event.particles_end(); ++particle) if (function(**particle)) return;
    print("no neutrino found");
}

Pythia8::Particle retrive_neutrino(HepMC::GenEvent const& event, double lifetime) {
    Pythia8::Particle pythia_particle;
    for_each_until(event, [&pythia_particle](HepMC::GenParticle const & hep_particle) {
        if (!is_heavy_neutral_lepton(hep_particle.pdg_id())) return false;
        auto const& momentum = hep_particle.momentum();
        pythia_particle = Pythia8::Particle(hep_particle.pdg_id(), hep_particle.status(), 0, 0, 0, 0, 0, 0, to_pythia(momentum), momentum.m(), momentum.m(), 9.);
        pythia_particle.vProd(to_pythia(hep_particle.production_vertex()->position()));
        return true;
    });
    pythia_particle.tau(lifetime);
    return pythia_particle;
}

void for_each_until(HepMC::IO_GenEvent& events, std::function<bool(HepMC::GenEvent const&)> const& function) {
    auto * event = events.read_next_event();
    if (!event) print("Hepmc file is empty");
    while (event) {
        if (function(*event)) break;
        delete event;
        events >> event;
    }
}

void for_each_until(Pythia8::Event const& event, std::function<void(Pythia8::Particle const& particle)> const& function) {
    for (auto line = 0; line < event.size(); ++line) function(event[line]);
}

double read(boost::filesystem::path const& path, Meta const& meta, double coupling) {
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
            if (!particle.isFinal() || !particle.isCharged()) continue;
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

void read(boost::filesystem::path const& path, double coupling) {
    auto meta = meta_info(path);
    if(!meta) return;
    Result result;
    result[meta->mass][coupling] = read(path, *meta, coupling);
    save(result, path.stem().string() + "-" + std::to_string(coupling));
}

boost::optional<Result> scan_file(boost::filesystem::path const& path) {
    Result result;
    print("file", path);
    if (auto meta = meta_info(path)) for (auto coupling : log_range(1e-8, 1, 8)) result[meta->mass][coupling] = read(path, *meta, coupling);
    else return boost::none;
    return result;
}

void scan(boost::filesystem::path const& path) {
    if(auto result = scan_file(path)) save(*result, path.stem().string());
}

void scans(std::string const& path_name) {
    if (debug) print("read hep mcs in", path_name);
    Result results;
    for (auto const& file : files(path_name)) if (file.path().extension().string() == ".hep") if (auto result = scan_file(file.path())) results += *result;
    save(results);
}

}

}

