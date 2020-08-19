#pragma once

#include "read-file.hh"
#include "pythia.hh"
#include "geometry.hh"
#include "ResonanceWidths.hh"

namespace hnl {

auto to_cgal(Pythia8::Particle const& particle) -> cgal::Point {
    return {particle.xProd() / 1000, particle.yProd() / 1000, particle.zProd() / 1000}; // convert from mm to m
}

auto max(std::map<int, std::map<int, double>> const& couplings) {
    double max = 0.;
    for (auto const& inner : couplings) for (auto const& pair : inner.second) if (pair.second > max) max = pair.second;
    return max;
}

inline double tau_to_Gamma(double tau) { //mm/c->GeV
    return 1.97327E-13 / tau;
}

auto for_each(Pythia8::ParticleDataEntry const& particle, std::function<void(Pythia8::DecayChannel const&)> const& function) {
    for (auto channel_number : irange(particle.sizeChannels())) function(particle.channel(channel_number));
}

std::vector<std::string> decay_table(std::function<double (int id_heavy, int id_light)> const& coupling, double mass, int id) {
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_passive(pythia);
    set_pythia_init(pythia);
    pythia.setResonancePtr(new NeutrinoResonance(pythia, coupling, mass, id));
    pythia.init();
    auto const& particle = *pythia.particleData.particleDataEntryPtr(id);
    std::vector<std::string> result({std::to_string(id) + ":new = " + pythia_hnl_name(id) + " void " + std::to_string(particle.spinType()) + " " + std::to_string(particle.chargeType()) + " " + std::to_string(particle.colType()) + " " + std::to_string(mass) + " " + to_string(particle.mWidth()) + " " + to_string(mass / 2) + " " + std::to_string(mass * 2) + " " + to_string(particle.tau0())});
    for_each(particle, [&result, id](Pythia8::DecayChannel const & channel) {
        int onMode = 1;
        int meMode = 101;
        result.emplace_back(std::to_string(id) + ":addChannel = " + std::to_string(onMode) + " " + std::to_string(channel.bRatio()) + " " + std::to_string(meMode) + " " + std::to_string(channel.product(0)) + " " + std::to_string(channel.product(1)) + " " + std::to_string(channel.product(2)));
    });
    return result;
}

double read_lhe(boost::filesystem::path const& path, Meta const& meta, double coupling) {
    print(path.string(), "with", coupling);
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_read_lhe(pythia, path.string());

    auto couplings = [&meta, coupling](int heavy, int light) {
        return meta.couplings.at(heavy).at(light) > 0 ? coupling : 0.;
    };
    for (auto const& line : decay_table(couplings, meta.mass, heavy_neutrino)) pythia.readString(line);

    pythia.init();

    int total = 0;
    int good = 0;
    int events_max = 1e6;
    auto analysis = mapp::analysis();
    for (auto event_number = 0; event_number < events_max; ++event_number) {
        ++total;
        if (!pythia.next()) {
            if (pythia.info.atEndOfFile()) break;
            print("Pythia encountered a problem");
            continue;
        }
        for (auto line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            if (debug) if (is_heavy_neutral_lepton(particle.id())) print(particle.id(), particle.m0(), particle.tau(), particle.isFinal(), particle.vDec().pAbs(), particle.mWidth());
            if (!particle.isFinal() || !particle.isCharged()) continue;
            if (!analysis.is_inside(to_cgal(particle))) continue;
            if (debug) print("Hooray!", particle.vProd().pAbs());
            ++good;
            break;
        }
    }
    pythia.stat();
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

double read_lhe(boost::filesystem::path const& path, double coupling) {
    auto meta = meta_info_lhe(path);
    return meta ? read_lhe(path, *meta, coupling) : 0.;
}

ScanResult scan_lhe(boost::filesystem::path const& path) {
    ScanResult result;
    print("file", path);
    auto meta = meta_info_lhe(path);
    if (meta) for (auto coupling : log_range(1e-8, 1, 80)) result[meta->mass][coupling] = read_lhe(path, *meta, coupling);
    return result;
}

void scan_lhe(std::string const& path_name) {
    boost::filesystem::path path("./" + path_name);
    save_result(scan_lhe(path), path.stem().string());
}

auto get_range(boost::filesystem::path const& path) {
    return boost::make_iterator_range(boost::filesystem::directory_iterator(path), {});
}

void scan_lhes(std::string const& path_name) {
    print("read lhe in", path_name);
    ScanResult result;
    for (auto const& folder : get_range(path_name)) if (boost::filesystem::is_directory(folder)) for (auto const& file : get_range(folder.path())) if (file.path().extension().string() == ".lhe" || file.path().extension().string() == ".gz") result += scan_lhe(file.path());
    save_result(result);
}

}
