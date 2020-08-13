#pragma once

#include <boost/range/algorithm/max_element.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

#include "generic.hh"
#include "geometry.hh"
#include "ResonanceWidths.hh"
#include "Sigma.hh"

namespace hnl {

namespace {

const bool debug = false;

std::vector<int> heavy_neutral_leptons() {
    return {9900012, 9900014, 9900016};
}

std::vector<int> light_neutrinos() {
    return {12, 14, 16};
}


const int heavy_neutrino = 9900012;

bool is_heavy_neutral_lepton(int id) {
    return id == 9900012 || id == 9900014 || id == 9900016;
}

template<typename First, typename Second>
std::map<First, Second>& operator+=(std::map<First, Second>& left, std::map<First, Second> const& right) {
    for (auto const& pair : right) left[pair.first] += pair.second;
    return left;
}

}

auto lin_scale(double min, double max, int step, int steps) {
    return min + (max - min) * step / steps;
}

auto log_scale(double min, double max, int step, int steps) {
    return std::pow(10, lin_scale(std::log10(min), std::log10(max), step, steps));
}

auto log_range(double min, double max, int steps) {
    return transform(irange(steps + 1), [&](auto step) {
        return log_scale(min, max, step, steps);
    });
}

auto neutrino_coupling = [](double factor) {
    return [factor](int id_heavy, int id_light) -> double {
        if (id_light != 12 && id_light != 14 && id_light != 16) {
            print(id_light, "is not a light neutrino");
            return 0;
        }
        if (id_heavy != 9900012 && id_heavy != 9900014 && id_heavy != 9900016) {
            print(id_heavy, "is not a heavy neutrino");
            return 0;
        }
        auto id = id_heavy - 9900000 - id_light;
        return id == 0 ? factor : 0;
        return id == 0 || std::abs(id) == 2 || std::abs(id) == 4 ? factor : 0;
    };
};

void set_pythia_production(Pythia8::Pythia& pythia) {
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000.");
    pythia.readString("ResonanceWidths:minWidth = 1E-30");
}

void set_pythia_init(Pythia8::Pythia& pythia) {
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showMultipartonInteractions = off");
}

void set_pythia_passive(Pythia8::Pythia& pythia) {
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("ResonanceWidths:minWidth = 1E-30");
}

void set_pythia_next(Pythia8::Pythia& pythia) {
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowInfo = 0");
}

void set_pythia_branching_fractions(Pythia8::Pythia& pythia) {
    set_pythia_production(pythia);
    set_pythia_init(pythia);
    set_pythia_passive(pythia);
}

void set_pythia_sigma(Pythia8::Pythia& pythia) {
    set_pythia_production(pythia);
    set_pythia_init(pythia);
    set_pythia_next(pythia);
}

auto has_neutrino = [](auto const& channel) {
    for (auto heavy : heavy_neutral_leptons()) {
        if (channel.product(0) == heavy || channel.product(1) == heavy || channel.product(2) == heavy || channel.product(3) == heavy || channel.product(4) == heavy) return true;
    }
    return false;
};

struct Loop {
    Loop(double min, int steps_) : m_min(min), steps(steps_) {}
    double m_min;
    int steps;
    double mass(double max, int step) const {
        return log_scale(m_min, max, step, steps);
    }
};

using Result = std::map<int, std::map<std::tuple<int, int, int, int, int>, std::map<int, double>>>;

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

Result branching_fraction(Loop const& loop, double& mass_max, int source, int step) {
    print(step, "of", loop.steps);
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_branching_fractions(pythia);
    if (!is_heavy_neutral_lepton(source)) mass_max = pythia.particleData.m0(source);
    is_heavy_neutral_lepton(source) ? pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling(1), loop.mass(mass_max, step), source)) : pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling(1), source));
    pythia.particleData.particleDataEntryPtr(source)->rescaleBR();
    pythia.init();
    Result result;
    auto const& particle = *pythia.particleData.particleDataEntryPtr(source);
    for (auto pos = 0; pos < particle.sizeChannels(); ++pos) {
        auto channel = particle.channel(pos);
        auto ratio = channel.bRatio();
        if (ratio > 0. && (has_neutrino(channel) || is_heavy_neutral_lepton(source)))
            result[source][ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}][step] = ratio;
    }
    return result;
}

void write_branching_fractions(int source) {
    Result result;
    Loop loop(.1, 5);
    double mass_max = 5;
    for (auto step = 0; step <= loop.steps; ++step) result += branching_fraction(loop, mass_max, source, step);
    save_data(result, loop, mass_max, source);
}

void write_branching_fractions() {

//     std::vector<int> sources{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541, 443, 553};
//     std::vector<int> sources{431, 411, 421};
//     std::vector<int> sources{511, 521, 531, 541};
//     std::vector<int> sources{511, 521, 531};
//     std::vector<int> sources{443, 553};
    std::vector<int> sources{heavy_neutrino};
    for (auto source : sources) write_branching_fractions(source);
}

void set_pythia_write_hepmc(Pythia8::Pythia& pythia, double mass) {
    pythia.particleData.m0(heavy_neutrino, mass);
    pythia.particleData.mMin(heavy_neutrino, 0.);
    pythia.particleData.particleDataEntryPtr(heavy_neutrino)->clearChannels();
    set_pythia_production(pythia);
    set_pythia_next(pythia);
    pythia.readString("Bottomonium:all = on");
    if (mass < 1.96849) pythia.readString("Charmonium:all = on");
    if (mass < .493677) pythia.readString("HardQCD:all = on");
//     pythia.readString("Onia:all(3S1) = on");
//     pythia.readString("SoftQCD:nonDiffractive = on");
//     pythia.readString("HardQCD:qqbar2ccbar  = on");
//     pythia.readString("HardQCD:gg2bbbar  = on");
//     pythia.readString("PhaseSpace:pTHatMin = .1");
    pythia.readString("Main:numberOfEvents = 100000");
}

std::string to_string(double value) {
    std::stringstream string_stream;
    string_stream << std::scientific << value;
    return string_stream.str();
}

void for_each_if(Pythia8::Pythia& pythia, HepMC::IO_GenEvent& hepmc_file, std::function<bool(Pythia8::Event const& event)> const& function) {
    HepMC::Pythia8ToHepMC converter;
    int event_number = 0;
    while (event_number < pythia.mode("Main:numberOfEvents")) {
        if (!pythia.next()) continue;
        if (!function(pythia.event)) continue;
        ++event_number;
        HepMC::GenEvent gen_event;
        converter.fill_next_event(pythia, &gen_event);
        hepmc_file.write_event(&gen_event);
    }
}

std::pair<int, double> get_max_width(Pythia8::Pythia& pythia, std::vector<int> const& mesons) {
    double max_width = 0;
    int id;
    for (auto meson : mesons) {
        double partial_width = 0.;
        auto& particle = *pythia.particleData.particleDataEntryPtr(meson);
        for (auto number : irange(particle.sizeChannels())) {
            auto& channel = particle.channel(number);
            if (has_neutrino(channel) && channel.bRatio() > 0.) {
                partial_width += channel.bRatio();
                if (debug) print(meson, channel.product(0), channel.product(1), channel.product(2));
            }
        }
        if (max_width < partial_width) {
            max_width = partial_width;
            id = meson;
        }
    }
    return {id, max_width};
}

auto get_optimal_coupling(double mass, std::vector<int> const& mesons) {
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_write_hepmc(pythia, mass);
    set_pythia_init(pythia);
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

void write_hepmc(double mass) {
    print("Generating events for HNLs with mass", mass, "GeV");

    std::vector<int> mesons{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541, 443, 553};
    auto coupling = get_optimal_coupling(mass, mesons);

    print("Use U^2 =", coupling);

    HepMC::IO_GenEvent hepmc_file(std::to_string(mass) + ".hep", std::ios::out);
    hepmc_file.write_comment("mass " + std::to_string(mass) + " GeV");


    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_write_hepmc(pythia, mass);
    for (auto meson : mesons) pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling(coupling), meson));

    pythia.init();

    auto [id, max_width] = get_max_width(pythia, mesons);
    print("The maximal BR into HNLs with U^2 =", coupling, "is", max_width, "from", id);

    for (auto heavy : heavy_neutral_leptons()) for (auto light : light_neutrinos()) hepmc_file.write_comment("coupling " + std::to_string(heavy) + " " + std::to_string(light) + " " + std::to_string(neutrino_coupling(coupling)(heavy, light)));

    int total = 0;
    int successfull = 0;
    int too_many = 0;

    for_each_if(pythia, hepmc_file, [&too_many, &total, &successfull](Pythia8::Event const & event) -> bool {
        ++total;
        int found_one = 0;
        bool success = false;
        for (auto line = 0; line < event.size(); ++line) {
            auto const& particle = event[line];
            if (!is_heavy_neutral_lepton(std::abs(particle.id()))) continue;
            ++found_one;
            if (found_one > 1) {
                print("More than one neutrino in event", total, "in line", line, "from mother", event[particle.mother1()].id());
                success = false;
                break;
            }
            double eta = std::abs(particle.eta());
            if (!(1.4 < eta && eta < 3.5)) continue;
            print("Success event", successfull + 1, "of", total, "in line", line, "from mother", event[particle.mother1()].id(), "with decay vertex", particle.vDec().pAbs(), "and eta", particle.eta(), "and phi", particle.phi());
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

    hepmc_file.write_comment("sigma " + to_string(pythia.info.sigmaGen() * successfull / total) + " mb");
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

Pythia8::Vec4 to_pythia(HepMC::FourVector const& vector) {
    return {vector.x(), vector.y(), vector.z(), vector.t()};
}

void for_each_until(HepMC::GenEvent const& gen_event, std::function<bool(HepMC::GenParticle const&)> const& function) {
    for (auto iterator = gen_event.particles_begin(); iterator != gen_event.particles_end(); ++iterator) if (function(**iterator)) return;
    print("no neutrino found");
}

auto retrive_neutrino(HepMC::GenEvent const& gen_event, double lifetime) -> Pythia8::Particle {
    Pythia8::Particle pythia_particle;
    for_each_until(gen_event, [&pythia_particle, &lifetime](HepMC::GenParticle const & hep_particle) {
        if (!is_heavy_neutral_lepton(hep_particle.pdg_id())) return false;
        auto& momentum = hep_particle.momentum();
        pythia_particle = Pythia8::Particle(hep_particle.pdg_id(), hep_particle.status(), 0, 0, 0, 0, 0, 0, to_pythia(momentum), momentum.m(), momentum.m(), 9.);
        pythia_particle.vProd(to_pythia(hep_particle.production_vertex()->position()));
        pythia_particle.tau(lifetime);
        return true;
    });
    return pythia_particle;
}

struct Line {
    friend std::istream& operator>>(std::istream& stream, Line& line) noexcept {
        std::getline(stream, line.string);
        return stream;
    }
    operator std::string() const noexcept {
        return string;
    }
private:
    std::string string;
};

auto split_line(std::string const& line) noexcept {
    std::vector<std::string> strings;
    boost::split(strings, line, [](char c) noexcept {
        return c == ' ';
    }, boost::token_compress_on);
    return strings;
}

auto import_file(boost::filesystem::path const& path) noexcept {
    std::ifstream file(path.string());
    std::vector<std::string> lines;
    std::copy(std::istream_iterator<Line>(file), std::istream_iterator<Line>(), std::back_inserter(lines));
    return lines;
}

template<typename Predicate>
auto read_file(std::vector<std::string>& lines, int pos, Predicate predicate) noexcept {
    auto found = boost::range::find_if(lines, [&predicate](auto & line) noexcept {
        boost::trim_if(line, boost::is_any_of("\t "));
        return predicate(split_line(line));
    });
    return found == lines.end() ? "value not found"s : split_line(*found).at(pos);
}

template<typename Predicate>
auto read_file(boost::filesystem::path const& path, int pos, Predicate predicate) noexcept {
    std::ifstream file(path.string());
    std::vector<std::string> lines;
    std::copy(std::istream_iterator<Line>(file), std::istream_iterator<Line>(), std::back_inserter(lines));
    auto found = boost::range::find_if(lines, [&predicate](auto & line) noexcept {
        boost::trim_if(line, boost::is_any_of("\t "));
        return predicate(split_line(line));
    });
    return found == lines.end() ? "value not found"s : split_line(*found).at(pos);
}

double to_double(std::string const& string) {
    try {
        return std::stod(string);
    } catch (...) {
        print("The string:", string, ", is not a number");
        return 0.;
    }
}

auto find_sigma(std::vector<std::string>& path) {
    if (debug) print("find sigma");
    return read_file(path, 1, [](auto const & strings)  {
        return strings.size() == 3 && strings.at(0) == "sigma" && strings.at(2) == "mb";
    });
}

auto find_mass(std::vector<std::string>& path) {
    if (debug) print("find mass");
    return read_file(path, 1, [](auto const & strings)  {
        return strings.size() == 3 && strings.at(0) == "mass" && strings.at(2) == "GeV";
    });
}

auto find_coupling(std::vector<std::string>& path, int heavy, int light) {
    if (debug) print("find coupling");
    return read_file(path, 3, [&](auto const & strings)  {
        return strings.size() == 4 && strings.at(0) == "coupling" && strings.at(1) == std::to_string(heavy) && strings.at(2) == std::to_string(light);
    });
}

auto to_cgal(Pythia8::Particle const& particle) -> cgal::Point {
    return {particle.xProd() / 1000, particle.yProd() / 1000, particle.zProd() / 1000}; // convert from mm to m
}

void set_pythia_read_hepmc(Pythia8::Pythia& pythia) {
    set_pythia_init(pythia);
    set_pythia_next(pythia);
    set_pythia_passive(pythia);
}

void for_each_until(HepMC::IO_GenEvent& hepmc_file, std::function<bool(HepMC::GenEvent const&)> const& function) {
    auto* hepmc_event = hepmc_file.read_next_event();
    if (!hepmc_event) print("Hepmc file is empty");
    while (hepmc_event) {
        if (!function(*hepmc_event)) break;
        delete hepmc_event;
        hepmc_file >> hepmc_event;
    }
}

struct Meta {
    double mass = 0;
    double sigma = 0;
    std::map<int, std::map<int, double>> couplings;
};

auto max(std::map<int, std::map<int, double>> const& couplings) {
    double max = 0.;
    for (auto const& inner : couplings) for (auto const& pair : inner.second) if (pair.second > max) max = pair.second;
    return max;
}

double read_hepmc(boost::filesystem::path const& path, Meta const& meta, double coupling) {
    if (debug) print("read hep mc", path.string(), "with", coupling);
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_read_hepmc(pythia);

    auto coupling_function = [&meta, coupling](int heavy, int light) {
        return meta.couplings.at(heavy).at(light) > 0 ? coupling : 0.;
    };
    pythia.setResonancePtr(new NeutrinoResonance(pythia, coupling_function, meta.mass, heavy_neutrino));
    pythia.init();

    auto lifetime = pythia.particleData.tau0(heavy_neutrino);
    if (debug) print("trying to open", path.string());
    HepMC::IO_GenEvent hepmc_file(path.string(), std::ios::in);
    if (debug) print("with result", hepmc_file.error_message());
    int total = 0;
    int good = 0;
    auto analysis = mapp::analysis();
    for_each_until(hepmc_file, [&](HepMC::GenEvent const & hepmc_event) -> bool {
        ++total;
        pythia.event.reset();
        pythia.event.append(retrive_neutrino(hepmc_event, lifetime));
        if (!pythia.next()) {
            print("Pythia encountered a problem");
            return true;
        }
        if (debug) pythia.event.list(true);
        for (auto line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            auto vertex = particle.vProd();
            if (particle.chargeType() == 0) continue;
            if (!analysis.is_inside(to_cgal(particle))) continue;
            if (debug) print("Hooray!", vertex.pAbs());
            ++good;
            break;
        }
        return true;
        if (total > 100) return false;
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

boost::optional<Meta> meta_info(boost::filesystem::path const& path) {
    auto file = import_file(path);
    Meta meta;
    meta.mass = to_double(find_mass(file));
    if (meta.mass <= 0) return boost::none;
    meta.sigma = to_double(find_sigma(file));
    if (meta.sigma <= 0) return boost::none;
    for (auto heavy : heavy_neutral_leptons()) for (auto light : light_neutrinos()) meta.couplings[heavy][light] = to_double(find_coupling(file, heavy, light));
    if (meta.couplings.empty()) return boost::none;
    return meta;
}

double read_hepmc(boost::filesystem::path const& path, double coupling) {
    auto meta = meta_info(path);
    return meta ? read_hepmc(path, *meta, coupling) : 0.;
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

ScanResult scan_hepmc(boost::filesystem::path const& path) {
    ScanResult result;
    print("file", path);
    auto meta = meta_info(path);
    if (meta) for (auto coupling : log_range(1e-6, 1, 6)) result[meta->mass][coupling] = read_hepmc(path, *meta, coupling);
    return result;
}

void scan_hepmc(std::string const& path_name) {
    boost::filesystem::path path("./" + path_name);
    save_result(scan_hepmc(path), path.stem().string());
}

void scan_hepmcs(std::string const& path_name) {
    if (debug) print("read hep mcs in", path_name);
    ScanResult result;
    for (auto const& directory_entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path_name), {})) if (directory_entry.path().extension().string() == ".hep") result += scan_hepmc(directory_entry.path());
    save_result(result);
}

void calculate_sigma(double mass) {
    if (debug) print("calculate_sigma");
    auto coupling = 1;
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_sigma(pythia);
    pythia.readString("Main:numberOfEvents = 100");
    pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling(coupling), mass, heavy_neutrino));
    Pythia8::SigmaProcess* sigma = new Sigma(heavy_neutrino, coupling);
    pythia.setSigmaPtr(sigma);
    pythia.init();
    for (auto event_number = 0; event_number < pythia.mode("Main:numberOfEvents"); ++event_number) if (!pythia.next()) print("Error in event", event_number);
    pythia.stat();
    delete sigma;
    print("sigma", pythia.info.sigmaGen() / coupling);
}

void write_sigma_hepmc(double mass) {
    print("Generating events for HNLs with mass", mass, "GeV");
    auto coupling = 1;
    print("Use U^2 =", coupling);

    HepMC::IO_GenEvent hepmc_file(std::to_string(mass) + ".hep", std::ios::out);
    hepmc_file.write_comment("mass " + std::to_string(mass) + " GeV");

    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_sigma(pythia);
    pythia.readString("Main:numberOfEvents = 100000");

    pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling(coupling), mass, heavy_neutrino));
    Pythia8::SigmaProcess* sigma = new Sigma(heavy_neutrino, coupling);
    pythia.setSigmaPtr(sigma);

    pythia.init();

    for (auto heavy : heavy_neutral_leptons()) for (auto light : light_neutrinos()) hepmc_file.write_comment("coupling " + std::to_string(heavy) + " " + std::to_string(light) + " " + std::to_string(neutrino_coupling(coupling)(heavy, light)));

    int total = 0;
    int successfull = 0;
    int too_many = 0;

    for_each_if(pythia, hepmc_file, [&too_many, &total, &successfull](Pythia8::Event const & event) -> bool {
        ++total;
        int found_one = 0;
        bool success = false;
        for (auto line = 0; line < event.size(); ++line) {
            auto const& particle = event[line];
            if (!is_heavy_neutral_lepton(std::abs(particle.id()))) continue;
            if(!is_heavy_neutral_lepton(event[particle.mother1()].id())) ++found_one;
            if (found_one > 1) {
                print("More than one neutrino in event", total, "in line", line, "from mother", event[particle.mother1()].id());
                success = false;
                break;
            }
            double eta = std::abs(particle.eta());
            if (!(1.4 < eta && eta < 3.5)) continue;
            if(!is_heavy_neutral_lepton(event[particle.mother1()].id())) print("Success event", successfull + 1, "of", total, "in line", line, "from mother", event[particle.mother1()].id(), "with decay vertex", particle.vDec().pAbs(), "and eta", particle.eta(), "and phi", particle.phi());
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

    hepmc_file.write_comment("sigma " + to_string(pythia.info.sigmaGen() * successfull / total) + " mb");
    delete sigma;
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

}

