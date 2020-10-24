#pragma once

#include <algorithm>

#include "read-file.hh"
#include "pythia-cgal.hh"
#include "mapp.hh"
#include "range.hh"
#include "ResonanceWidths.hh"
#include "hepmc.hh"
#include "string.hh"

namespace hnl {

namespace hepmc {

namespace {

const bool debug = true;

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
    if (debug) print("meta info", path.string());
//     auto lines = import_head(path, 100) + import_tail(path, 100);
    auto a = import_head(path, 100);
    if (debug) print(a.size());
    if (debug) print(a);
    auto b = import_tail(path, 100);
    if (debug) print(b.size());
    if (debug) print(b);
    auto lines = a + b;
    if (debug) print(lines.size());
    if (debug) print(lines);
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

HepMC::GenParticle retrive_neutrino_2(HepMC::GenEvent const& event) {
    if (debug) print("get neutrino");
    HepMC::GenParticle pythia_particle;
    for_each_until(event, [&pythia_particle](HepMC::GenParticle const & hep_particle) {
        if (!is_heavy_neutral_lepton(hep_particle.pdg_id())) return false;
        pythia_particle = hep_particle;
        return true;
    });
    return pythia_particle;
}

void for_each_until(HepMC::IO_GenEvent& events, std::function<bool(HepMC::GenEvent const&)> const& function) {
    auto* event = events.read_next_event();
    if (!event) print("Hepmc file is empty");
    int counter = 0;
    while (event) {
        ++counter;
        if (debug) print("Event", counter);
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
    long long total = 0;
    long long good = 0;
    long long events_max = LLONG_MAX;
    long long max_tries = 10e5;
    auto analysis = mapp::analysis();
    for_each_until(events, [&](HepMC::GenEvent const & event) -> bool {
        auto neutrino = retrive_neutrino(event, lifetime);
        bool is_good = false;
        long long sub_total = 0;
        do {
            ++sub_total;
            pythia.event.reset();
            pythia.event.append(neutrino);
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
                is_good = true;
                break;
            }
        } while (sub_total < max_tries && !is_good);
        total += sub_total;
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
    if (!meta) return;
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
    if (auto result = scan_file(path)) save(*result, path.stem().string());
}

void scans(std::string const& path_name) {
    if (debug) print("read hep mcs in", path_name);
    Result results;
    for (auto const& file : files(path_name)) if (file.path().extension().string() == ".hep") if (auto result = scan_file(file.path())) results += *result;
    save(results);
}

double three(HepMC::FourVector const& vector) {
    return std::sqrt(sqr(vector.x()) + sqr(vector.y()) + sqr(vector.z()));
}

double betagamma(HepMC::FourVector const& vector) {
    return vector.m() / three(vector);
}

// auto russian() {
//     double l1 = 5; // m
//     double l2 = 7; // m
//     double beta = 1; // m
//     return std::exp(l1 / beta) - std::exp(l2 / beta);
// }

std::vector<std::pair<double, int>> histogram(std::vector<double> const& data) {
    if (debug) print("histogram");
    auto const [min, max] = std::minmax_element(begin(data), end(data));
    int bins = 100;
    std::vector<std::pair<double, int>> histogram(bins, {0., 0});
    for (auto i = 0; i < bins; ++i) histogram[i].first = log_value(*min, *max, i, bins);
    for (auto point : data) {
        int i = static_cast<int>(std::floor(log_step(*min, *max, point, bins)));
        if (i == bins) --i;
        if (i < 0 || i >= bins) print("going to acces", i);
        histogram[i].second++;
    }
    return histogram;
}

std::vector<std::pair<double, int>> read_simplified_det(boost::filesystem::path const& path, double coupling) {
    if (debug) print("read hep mc", path.string(), "with", coupling);

//     auto couplings = [&meta, coupling](int heavy, int light) {
//         return meta.couplings.at(heavy).at(light) > 0 ? coupling : 0.;
//     };

    if (debug) print("trying to open", path.string());
    HepMC::IO_GenEvent events(path.string(), std::ios::in);
    if (debug) print("with result", events.error_message());

    std::vector<double> betas;
    std::size_t events_max = 10000000;

    for_each_until(events, [&](HepMC::GenEvent const & event) -> bool {
        auto neutrino = retrive_neutrino_2(event);
        betas.emplace_back(betagamma(neutrino.momentum()));
        return betas.size() > events_max;
    });


    return histogram(betas);
}

void save(std::vector<std::pair<double, int>> const& result, std::string const& name) {
    std::ofstream file;
    file.open(name + ".dat");
    bool first = false;
    for (auto const& line : result) {
        if (first) {
            file << "mass" << '\t';
            file << std::scientific << line.first << '\t' << line.second;
            file << std::endl;
            first = false;
        }
        file << std::scientific << line.first << '\t' << line.second;
        file << std::endl;
    }
}

void read_simplified(boost::filesystem::path const& path, double coupling) {
    if (debug) print("read hep mc simp");
    auto result = read_simplified_det(path, coupling);
    save(result, path.stem().string() + "log");
}

void save(std::vector<Meta> const& metas, std::string const& name) {
    std::ofstream file;
    file.open(name + ".dat");
    bool first = true;
    for (auto const& meta : metas) {
        if (first) {
            file << "mass [GeV]" << '\t' << "crosssection [mb]" << '\t' << "coupling" << '\t' << "scaled sigma" << std::endl;
            first = false;
        }
        file << std::scientific << meta << '\t' << meta.sigma / max(meta.couplings) << std::endl;
    }
}

void extract_metas(boost::filesystem::path const& path) {
    if (debug) print("extract meta", path.string());
    std::vector<Meta> metas;
    for (auto const& file : files(path)) if (file.path().extension().string() == ".hep" || file.path().extension().string() == ".gz") if (auto meta = meta_info(file)) metas.emplace_back(*meta);
    save(metas, "meta");
}

}

}

