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


double read_lhe(boost::filesystem::path const& path, Meta const& meta, double coupling) {
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);

    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = " + path.string());

    auto couplings = [&meta, coupling](int heavy, int light) {
        return meta.couplings.at(heavy).at(light) > 0 ? coupling : 0.;
    };
    pythia.setResonancePtr(new NeutrinoResonance(pythia, couplings, meta.mass, heavy_neutrino));
    pythia.init();

//     auto lifetime = pythia.particleData.tau0(heavy_neutrino);

    int total = 0;
    int good = 0;
//     int events_max = 1e6;
    auto analysis = mapp::analysis();
    for (auto event_number = 0; ; ++event_number) {
        ++total;
        if (!pythia.next()) {
            if(pythia.info.atEndOfFile()) break;
            print("Pythia encountered a problem");
            continue;
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
    if (meta) for (auto coupling : log_range(1e-8, 1, 8)) result[meta->mass][coupling] = read_lhe(path, *meta, coupling);
    return result;
}

void scan_lhe(std::string const& path_name) {
    boost::filesystem::path path("./" + path_name);
    save_result(scan_lhe(path), path.stem().string());
}

auto get_range(boost::filesystem::path const& path){
    return boost::make_iterator_range(boost::filesystem::directory_iterator(path), {});
}

void scan_lhes(std::string const& path_name) {
    print("read lhe in", path_name);
    ScanResult result;
    for (auto const& folder : get_range(path_name)) if(boost::filesystem::is_directory(folder)) for (auto const& file : get_range(folder.path())) if (file.path().extension().string() == ".lhe" || file.path().extension().string() == ".gz") result += scan_lhe(file.path());
    save_result(result);
}

}
