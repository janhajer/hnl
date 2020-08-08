#include <boost/range/algorithm/max_element.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

#include "generic.hh"
#include "geometry.hh"
#include "ResonanceWidths.hh"

namespace {

const bool debug = true;

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

}

namespace hnl {

auto lin_scale(double min, double max, int step, int steps) {
    return min + (max - min) * step / steps;
}

auto log_scale(double min, double max, int step, int steps) {
    return std::pow(10, lin_scale(std::log10(min), std::log10(max), step, steps));
}

auto log_range(double min, double max, int steps) {
    return transform(irange(steps + 1), [&](auto step)  {
        return log_scale(min, max, step, steps);
    });
}

auto neutrino_coupling = [](int id_heavy, int id_light) -> double {
    if (id_light != 12 && id_light != 14 && id_light != 16) {
        print(id_light, "is not a light neutrino");
        return 0;
    }
    if (id_heavy != 9900012 && id_heavy != 9900014 && id_heavy != 9900016) {
        print(id_heavy, "is not a heavy neutrino");
        return 0;
    }
    auto id = id_heavy - 9900000 - id_light;
    return id == 0 ? 1e-1 : 0;
    return id == 0 || std::abs(id) == 2 || std::abs(id) == 4 ? .1 : 0;
};

void set_pythia_production(Pythia8::Pythia& pythia) {
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000.");
    pythia.readString("ResonanceWidths:minWidth = 1E-30");
//     pythia.particleData.m0(1, 0.0048);
//     pythia.particleData.m0(2, 0.0023);
//     pythia.particleData.m0(3, 0.095);
//     pythia.particleData.m0(4, 1.275);
//     pythia.particleData.m0(5, 4.180);
}

void set_pythia_init(Pythia8::Pythia& pythia) {
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Init:showProcesses = off");
    pythia.readString("Init:showChangedSettings = off");
}

void set_pythia_passive(Pythia8::Pythia& pythia) {
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("ResonanceWidths:minWidth = 1E-30");
}

void set_pythia_branching_fractions(Pythia8::Pythia& pythia) {
    set_pythia_production(pythia);
    set_pythia_init(pythia);
    set_pythia_passive(pythia);
}

auto has_neutrino = [](auto const& channel) {
    return channel.product(0) == heavy_neutrino || channel.product(1) == heavy_neutrino || channel.product(2) == heavy_neutrino || channel.product(3) == heavy_neutrino || channel.product(4) == heavy_neutrino;
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

int write_branching_fractions() {

//     std::vector<int> sources{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541, 443, 553};
//     std::vector<int> sources{431, 411, 421};
//     std::vector<int> sources{511, 521, 531, 541};
//     std::vector<int> sources{511, 521, 531};
//     std::vector<int> sources{443, 553};
    std::vector<int> sources{heavy_neutrino};
    Result result;
    for (auto source : sources) {
        Loop loop(.1, 5);
        double mass_max = 5;
        for (auto step = 0; step <= loop.steps; ++step) {
            print(step, "of", loop.steps);
            Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
            set_pythia_branching_fractions(pythia);
            if (!is_heavy_neutral_lepton(source)) mass_max = pythia.particleData.m0(source);
            pythia.particleData.m0(heavy_neutrino, loop.mass(mass_max, step));
            is_heavy_neutral_lepton(source) ? pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling, source)) : pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling, source));
            pythia.particleData.particleDataEntryPtr(source)->rescaleBR();
            pythia.init();
            auto const& particle = *pythia.particleData.particleDataEntryPtr(source);
            for (auto pos = 0; pos < particle.sizeChannels(); ++pos) {
                auto channel = particle.channel(pos);
                auto ratio = channel.bRatio();
                if (ratio > 0. && (has_neutrino(channel) || is_heavy_neutral_lepton(source)))
                    result[source][ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}][step] = ratio;
            }
        }
        save_data(result, loop, mass_max, source);
    }
    return 0;
}

void set_pythia_next(Pythia8::Pythia& pythia) {
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberShowProcess = 0");
//     pythia.readString("Next:numberShowInfo = 0");
}

void set_pythia_write_hepmc(Pythia8::Pythia& pythia, double mass) {
    pythia.particleData.m0(heavy_neutrino, mass);
    pythia.particleData.mMin(heavy_neutrino, 0.);
    set_pythia_production(pythia);
    set_pythia_next(pythia);
    pythia.readString("Bottomonium:all = on");
    if (mass < 1.86962) pythia.readString("Charmonium:all = on");
    if (mass < .493677) pythia.readString("HardQCD:all = on");
//     pythia.readString("Onia:all(3S1) = on");
//     pythia.readString("SoftQCD:nonDiffractive = on");
//     pythia.readString("HardQCD:qqbar2ccbar  = on");
//     pythia.readString("HardQCD:gg2bbbar  = on");
//     pythia.readString("PhaseSpace:pTHatMin = .1");
    pythia.readString("Main:numberOfEvents = 100000");
}

int write_hepmc(double mass) {
    HepMC::IO_GenEvent hepmc_file("neutrino_" + std::to_string(mass) + ".hep", std::ios::out);
    hepmc_file.write_comment("mass " + std::to_string(mass) + " GeV");

    for (auto heavy : heavy_neutral_leptons()) for (auto light : light_neutrinos()) hepmc_file.write_comment("Coupling " + std::to_string(heavy) + " " + std::to_string(light) + " " + std::to_string(neutrino_coupling(heavy, light)) + " GeV");

    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_write_hepmc(pythia, mass);
    std::vector<int> mesons{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541};
    for (auto meson : mesons) pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling, meson));
    pythia.init();

    HepMC::Pythia8ToHepMC converter;
    int too_many = 0;
    int total = 0;
    int successfull = 0;
    while (successfull < pythia.mode("Main:numberOfEvents")) {
        if (!pythia.next()) continue;
        ++total;
        int found_one = 0;
        bool success = false;
        for (auto line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            if (!is_heavy_neutral_lepton(std::abs(particle.id()))) continue;
            ++found_one;
            if (found_one > 1) {
                print("More than one neutrino in event", total);
                success = false;
                break;
            }
            double eta = std::abs(particle.eta());
            if (!(1.4 < eta && eta < 3.5)) continue;
            print("Success event", successfull, "of", total, "in line", line, "from mother", pythia.event[particle.mother1()].id(), "with decay vertex", particle.vDec().pAbs(), "and eta", particle.eta(), "and phi", particle.phi());
            success = true;
        }
        if (found_one > 1) ++too_many;
        if (!success) continue;
        ++successfull;
        HepMC::GenEvent gen_event;
        converter.fill_next_event(pythia, &gen_event);
        hepmc_file.write_event(&gen_event);
    }

    pythia.stat();

    print("\nwritten", successfull, "events of", total, "that is a fraction of", double(successfull) / total);
    print("\ntherefore from", pythia.info.sigmaGen(), "mb we save", pythia.info.sigmaGen() * successfull / total, "mb");
    print("\nfound", too_many, "events with more than one neutrino. that is a fraction of", double(too_many) / total, "and", double(too_many) / successfull);

    hepmc_file.write_comment("sigma " + std::to_string(pythia.info.sigmaGen() * successfull / total) + " mb");
    return 0;
}

int write_hepmcs() {
    Loop loop(.1, 50);
    for (auto step = 0; step <= loop.steps; ++step) {
        auto mass = loop.mass(6., step);
        std::ofstream ofstream("neutrino_" + std::to_string(mass) + ".txt");
        std::streambuf* streambuf = std::cout.rdbuf();
        std::cout.rdbuf(ofstream.rdbuf());
        write_hepmc(mass);
        std::cout.rdbuf(streambuf);
    }
    return 0;
}

auto retrive_neutrino(HepMC::GenEvent const* const gen_event, double lifetime) -> Pythia8::Particle {
    for (auto iterator = gen_event->particles_begin(); iterator != gen_event->particles_end(); ++iterator) {
        auto const& particle = **iterator;
        if (!is_heavy_neutral_lepton(particle.pdg_id())) continue;
        Pythia8::Particle part(particle.pdg_id(), particle.status(), 0, 0, 0, 0, 0, 0, particle.momentum().px(), particle.momentum().py(), particle.momentum().pz(), particle.momentum().e(), particle.momentum().m(), particle.momentum().m(), 9.);
        part.xProd(particle.production_vertex()->position().x());
        part.yProd(particle.production_vertex()->position().y());
        part.zProd(particle.production_vertex()->position().z());
        part.tProd(particle.production_vertex()->position().t());
        part.tau(lifetime);
        return part;
    }
    print("no neutrino found");
    return {};
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

double convert(std::string const& string) {
    try {
        return std::stod(string);
    } catch (...) {
        print(string, "is not a number");
        return 0.;
    }
}

auto find_sigma(boost::filesystem::path const& path) {
    return read_file(path, 1, [](auto const & strings)  {
        return strings.size() == 3 && strings.at(0) == "sigma" && strings.at(2) == "mb";
    });
}

auto find_mass(boost::filesystem::path const& path) {
    return read_file(path, 1, [](auto const & strings)  {
        return strings.size() == 3 && strings.at(0) == "mass" && strings.at(2) == "GeV";
    });
}

auto find_coupling(boost::filesystem::path const& path, int heavy, int light) {
    return read_file(path, 3, [&](auto const & strings)  {
        return strings.size() == 4 && strings.at(0) == "coupling" && strings.at(1) == std::to_string(heavy) && strings.at(2) == std::to_string(light);
    });
}

auto get_point(Pythia8::Particle const& particle) -> cgal::Point {
    return {particle.xProd() / 1000, particle.yProd() / 1000, particle.zProd() / 1000}; // convert from mm to m
}

void set_pythia_read_hepmc(Pythia8::Pythia& pythia) {
    set_pythia_init(pythia);
    set_pythia_next(pythia);
    set_pythia_passive(pythia);
}

void for_each(HepMC::IO_GenEvent& hepmc_file, std::function<bool(HepMC::GenEvent const* const)> const& function) {
    auto* hepmc_event = hepmc_file.read_next_event();
    while (hepmc_event) {
        if (!function(hepmc_event)) break;
        delete hepmc_event;
        hepmc_file >> hepmc_event;
    }
}

auto max(std::map<int, std::map<int, double>> const& couplings) {
    double max = 0.;
    for (auto const& inner : couplings) for (auto const& pair : inner.second) if (pair.second > max) max = pair.second;
    return max;
}

double read_hepmc(boost::filesystem::path path, double factor = 1.) {
    if(debug) print("read hep mc", path.string(), "with", factor);
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    set_pythia_read_hepmc(pythia);
    print("get values");
    auto mass = convert(find_mass(path));
    auto sigma = convert(find_sigma(path));
    print("got values");
    if (sigma <= 0.) return 0.;
    pythia.particleData.m0(heavy_neutrino, mass);

    std::map<int, std::map<int, double>> couplings;
    for (auto heavy : heavy_neutral_leptons()) for (auto light : light_neutrinos()) couplings[heavy][light] = convert(find_coupling(path, heavy, light));
    auto coupling_function = [&couplings, factor](int heavy, int light) {
        return couplings[heavy][light] * factor;
    };
    pythia.setResonancePtr(new NeutrinoResonance(pythia, coupling_function, heavy_neutrino));
    pythia.init();

    auto lifetime = pythia.particleData.tau0(heavy_neutrino);
    HepMC::IO_GenEvent hepmc_file(path.filename().string(), std::ios::in);
    int total = 0;
    int good = 0;
    auto analysis = mapp::analysis();
    print("lets go");
    for_each(hepmc_file, [&](HepMC::GenEvent const * const hepmc_event) -> bool {
        ++total;
        pythia.event.reset();
        pythia.event.append(retrive_neutrino(hepmc_event, lifetime));
        if (!pythia.next()) {
            print("Pythia encountered a problem");
            return true;
        }
        if (debug) pythia.event.list(true);

        print("total",total);
        if (total > 10) return false;

        for (auto line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            auto vertex = particle.vProd();
            if (particle.chargeType() == 0) continue;
            if (!analysis.is_inside(get_point(particle))) continue;
            if (debug) print("Hooray!", vertex.pAbs());
            ++good;
            break;
        }
        return true;
    });
    print(total);
    auto fraction = double(good) / total;
    print(good, total, fraction);
    print("mass", mass, "GeV", "factor", factor, "coupling", max(couplings) * factor, "sigma", sigma * fraction * factor, "mb");
    return sigma * fraction * factor;
}

void save_result(std::map<double, std::map<double, double>> const& result) {
    std::ofstream file;
    file.open("result.dat");
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

int read_hepmcs(std::string const& path) {
    if(debug) print("read hep mcs", path);
    std::map<double, std::map<double, double>> result;
    for (auto const& file : boost::make_iterator_range(boost::filesystem::directory_iterator(path), {})) {
        if (file.path().extension().string() != ".hep") continue;
        auto mass = convert(find_mass(file.path()));
        if (mass <= 0) continue;
        std::map<int, std::map<int, double>> couplings;
        for (auto heavy : heavy_neutral_leptons()) for (auto light : light_neutrinos()) {
                auto value = convert(find_coupling(file.path(), heavy, light));
                if (value > 0) couplings[heavy][light] = value;
            }
        for (auto factor : log_range(1e-6, 1, 6)) result[mass][max(couplings) * factor] = read_hepmc(file.path(), factor);
    }
    save_result(result);
    return 0;
}

}

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    using namespace hnl;
    return read_hepmcs(arguments.empty() ? "." : arguments.front());
    return write_hepmc(arguments.empty() ? 1. : convert(arguments.front()));
    return write_branching_fractions();
    return read_hepmc(arguments.empty() ? "neutrino_0.500000.hep" : arguments.front());
    return write_hepmcs();
}


