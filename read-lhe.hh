#pragma once

#include "decay_table.hh"
#include "read-file.hh"
#include "pythia-cgal.hh"
#include "mapp.hh"
#include "string.hh"

namespace hnl {

namespace lhe {

namespace {

const bool debug = false;

}

auto mass(std::vector<std::string> const& lines) {
     if(debug) print("mass");
    return find_if(lines, 1, [](auto const & strings) {
        return strings.size() > 2 && strings.at(0) == std::to_string(heavy_neutrino) && strings.at(2) == "#" && strings.at(3) == "mn1";
    });
}

auto sigma(std::vector<std::string> const& lines) {
     if(debug) print("sigma");
    return find_if(lines, 5, [](auto const & strings) {
        return strings.size() > 4 && strings.at(0) == "#" && strings.at(1) == "Integrated" && strings.at(2) == "weight" && strings.at(3) == "(pb)" && strings.at(4) == ":";
    });
}

std::string parameter_heavy(int heavy){
    switch (heavy){
        case 9900012 : return "n1";
        case 9900014 : return "n2";
        case 9900016 : return "n3";
        default : print("not a heavy neutrino");
    }
    return "";
}

std::string parameter_light(int light){
    switch (light){
        case 12 : return "e";
        case 14 : return "mu";
        case 16 : return "ta";
        default : print("not a heavy neutrino");
    }
    return "";
}

std::string parameter(int heavy, int light){
    return "v" + parameter_light(light) + parameter_heavy(heavy);
}

std::string coupling(std::vector<std::string> const& lines, int heavy, int light, int pos) {
     if(debug) print("coupling", heavy, light, pos);
    std::string name = parameter(heavy, light);
    return find_if(lines, 1, [&name, pos](auto const & strings) {
        return strings.size() > 3 && strings.at(0) == std::to_string(pos) && strings.at(2) == "#" && strings.at(3) == name;
    });
}

boost::optional<Meta> meta_info(boost::filesystem::path const& path) {
    auto lines = import_head(path, 500);
    Meta meta;
    meta.mass = to_double(mass(lines));
    if (meta.mass <= 0) return boost::none;
    meta.sigma = to_double(sigma(lines)) * 1e-12 / 1e-3 ; // from picobarn (MadGraph) to millibarn (Pythia)
    if (meta.sigma <= 0) return boost::none;
    int pos = 0;
    for (auto light : light_neutrinos()) for (auto heavy : heavy_neutral_leptons()) meta.couplings[heavy][light] = to_double(coupling(lines, heavy, light, ++pos));
    if (meta.couplings.empty()) return boost::none;
    print("Meta info",meta.mass, meta.sigma);
    return meta;
}

double read(boost::filesystem::path const& path, Meta const& meta, double coupling) {
    print(path.string(), "with", coupling);
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_read_lhe(pythia, path.string());

    auto couplings = [&meta, coupling](int heavy, int light) {
        return meta.couplings.at(heavy).at(light) > 0 ? coupling : 0.;
    };
    for (auto const& line : hnl_decay_table(couplings, meta.mass, heavy_neutrino)) {
        if(debug) print(line);
        pythia.readString(line);
    }

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

double read(boost::filesystem::path const& path, double coupling) {
    auto meta = meta_info(path);
    return meta ? read(path, *meta, coupling) : 0.;
}

Result scan(boost::filesystem::path const& path) {
    Result result;
    print("file", path);
    auto meta = meta_info(path);
    if (meta) for (auto coupling : log_range(1e-8, 1, 8)) result[meta->mass][coupling] = read(path, *meta, coupling);
    return result;
}

void scan(std::string const& path_name) {
    boost::filesystem::path path("./" + path_name);
    save(scan(path), path.parent_path().string());
//     save(scan(path), path.stem().string());
}

void scans(std::string const& path_name) {
    print("read lhe in", path_name);
    Result result;
    for (auto const& folder : files(path_name)) if (boost::filesystem::is_directory(folder)) for (auto const& file : files(folder.path())) if (file.path().extension().string() == ".lhe" || file.path().extension().string() == ".gz") result += scan(file.path());
    save(result);
}

}

}
