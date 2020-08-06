#include "generic.hh"
#include "ResonanceWidths.hh"
#include "read_files.hh"

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "geometry.hh"

namespace
{

const bool debug = false;

const int heavy_neutrino = 9900012;

bool is_neutrino(int id)
{
    return id == 9900012 || id == 9900014 || id == 9900016;
}

double tau_to_Gamma(double tau)//mm/c->GeV
{
    return 1.97327E-13 / tau;
}

}

namespace hnl
{

struct UserHook : public Pythia8::UserHooks {
    UserHook(double mass_) : mass(mass_) {}
    virtual bool canVetoPartonLevel() override
    {
        return true;
    }
    virtual bool doVetoPartonLevel(Pythia8::Event const& event) override
    {
        for (auto line = 0; line < event.size(); ++line) {
            auto const& particle = event[line];
            print(line, particle.name(), particle.e());
            if (particle.m0() > mass) {
                print("is heavy enough", particle.name());
                return false;
            }
        }
        print("Not found");
        return true;
    }
private:
    double mass;
};

auto log_scale(double min, double max, int bin, int bin_number) noexcept
{
    auto log_min = std::log10(min);
    auto log_max = std::log10(max);
    auto bin_width = (log_max - log_min) / bin_number;
    return std::pow(10, log_min + bin * bin_width);
}


auto lin_scale(double min, double max, int bin, int bin_number) noexcept
{
    auto bin_width = (max - min) / bin_number;
    return min + bin * bin_width;
}

void check_line()
{
    static int line = 0;
    using namespace hnl;
    print("line", line);
    ++line;
}

auto neutrino_coupling_2 = [](int id_heavy, int id_light) -> double
{
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

struct Loop {
    Loop(double min, int steps_) : m_min(min), steps(steps_) {}
    double m_min;
    int steps;
    double mass(double max, int step) const
    {
        return log_scale(m_min, max, step, steps);
    }
};

}

void set_pythia_global(Pythia8::Pythia& pythia)
{
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000.");
    pythia.particleData.m0(1, 0.0048);
    pythia.particleData.m0(2, 0.0023);
    pythia.particleData.m0(3, 0.095);
    pythia.particleData.m0(4, 1.275);
    pythia.particleData.m0(5, 4.180);
}

void set_pythia_single(Pythia8::Pythia& pythia, double mass)
{
    pythia.particleData.m0(heavy_neutrino, mass);
    pythia.particleData.mMin(heavy_neutrino, 0.);
    set_pythia_global(pythia);
    pythia.readString("Bottomonium:all = on");
    if (mass < 1.86962) pythia.readString("Charmonium:all = on");
    if (mass < .493677) pythia.readString("HardQCD:all = on");
//     pythia.readString("Onia:all(3S1) = on");
//     pythia.readString("SoftQCD:nonDiffractive = on");
//     pythia.readString("HardQCD:qqbar2ccbar  = on");
//     pythia.readString("HardQCD:gg2bbbar  = on");
//     pythia.readString("PhaseSpace:pTHatMin = .1");
    pythia.readString("Next:numberShowEvent = 0");
//     pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Main:numberOfEvents = 1000");
}

int write_hepmc(double mass)
{
    using namespace hnl;

    HepMC::Pythia8ToHepMC pythia_to_hep;
    pythia_to_hep.set_store_proc();
    pythia_to_hep.set_store_xsec();
    pythia_to_hep.set_store_pdf();
    HepMC::IO_GenEvent io_event("neutrino_" + std::to_string(mass) + ".hep", std::ios::out);

    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);

    set_pythia_single(pythia, mass);

    std::vector<int> mesons{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541};
    for (auto meson : mesons) pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling_2, meson));

    pythia.init();

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
//             if (particle.id() > 100) print(particle.name());
            if (!is_neutrino(std::abs(particle.id()))) continue;
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
        pythia_to_hep.fill_next_event(pythia, &gen_event);
        io_event.write_event(&gen_event);
    }

    pythia.stat();

    print("\nwritten", successfull, "events of", total, "that is a fraction of", double(successfull) / total);
    print("\ntherefore from", pythia.info.sigmaGen(), "mb we save", pythia.info.sigmaGen() * successfull / total, "mb");
    print("\nfound", too_many, "events with more than one neutrino. that is a fraction of", double(too_many) / total, "and", double(too_many) / successfull);

    io_event.write_comment("sigma " + std::to_string(pythia.info.sigmaGen() * successfull / total) + " mb");
    return 0;
}

int write_hepmcs()
{
    using namespace hnl;
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

void set_pythia(Pythia8::Pythia& pythia)
{
    set_pythia_global(pythia);
    pythia.readString("Bottomonium:all = on");
//     pythia.readString("Charmonium:all = on");

    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Init:showProcesses  = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythia.readString("PhaseSpace:showViolation = on");
    pythia.readString("PhaseSpace:increaseMaximum = on");
//     pythia.readString("PhaseSpace:showSearch = on");
}

auto has_neutrino = [](auto const& channel)
{
    return channel.product(0) == heavy_neutrino || channel.product(1) == heavy_neutrino || channel.product(2) == heavy_neutrino || channel.product(3) == heavy_neutrino || channel.product(4) == heavy_neutrino;
};

template<typename Data>
void save_data(Data& result, hnl::Loop const& loop, double mass, int source)
{
    using namespace hnl;
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

int write_branching_fractions()
{
    using Result = std::map<int, std::map<std::tuple<int, int, int, int, int>, std::map<int, double>>>;
    using namespace hnl;

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
            set_pythia(pythia);
            if (!is_neutrino(source)) mass_max = pythia.particleData.m0(source) * 1.01;
            pythia.particleData.m0(heavy_neutrino, loop.mass(mass_max, step));
            is_neutrino(source) ? pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling_2, source)) : pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling_2, source));
            pythia.particleData.particleDataEntryPtr(source)->rescaleBR();
            pythia.init();
            auto const& particle = *pythia.particleData.particleDataEntryPtr(source);
            for (auto pos = 0; pos < particle.sizeChannels(); ++pos) {
                auto channel = particle.channel(pos);
                auto ratio = channel.bRatio();
                if (ratio > 0. && (has_neutrino(channel) || is_neutrino(source)))
                    result[source][ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}][step] = ratio;
            }
        }
        save_data(result, loop, mass_max, source);
    }
    return 0;
}

using namespace hnl;

auto retrive_neutrino(HepMC::GenEvent* const gen_event, double tau) -> Pythia8::Particle
{
    for (auto iterator = gen_event->particles_begin(); iterator != gen_event->particles_end(); ++iterator) {
        auto const& particle = **iterator;
        if (!is_neutrino(particle.pdg_id())) continue;
        Pythia8::Particle part(particle.pdg_id(), particle.status(), 0, 0, 0, 0, 0, 0, particle.momentum().px(), particle.momentum().py(), particle.momentum().pz(), particle.momentum().e(), particle.momentum().m(), particle.momentum().m(), 9.);
        part.xProd(particle.production_vertex()->position().x());
        part.yProd(particle.production_vertex()->position().y());
        part.zProd(particle.production_vertex()->position().z());
        part.tProd(particle.production_vertex()->position().t());
        part.tau(tau);
        return part;
    }
    print("no neutrino found");
    return Pythia8::Particle();
}

auto find_sigma(boost::filesystem::path const& path) noexcept
{
    return read_file(path, 1, [](auto const & strings) noexcept {
        return strings.size() == 3 && strings.at(0) == "sigma" && strings.at(2) == "mb";
    });
}

auto get_point(Pythia8::Particle const& particle) -> cgal::Point
{
    return {particle.xProd() / 1000, particle.yProd() / 1000, particle.zProd() / 1000};
}

int read_hepmc(std::string const& file_name)
{
    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    pythia.readString("ProcessLevel:all = off"); // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
// Switch off automatic event listing in favour of manual.
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("ResonanceWidths:minWidth = 1E-30");
    pythia.particleData.m0(heavy_neutrino, .2);
    pythia.particleData.mayDecay(heavy_neutrino, true);
    pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling_2, heavy_neutrino));
    pythia.init();
    boost::filesystem::path path(file_name);
    auto sigma = std::stod(find_sigma(path));
    print(sigma);

    pythia.particleData.isResonance(heavy_neutrino, false);
    pythia.particleData.mayDecay(heavy_neutrino, true);
    auto tau = pythia.particleData.tau0(heavy_neutrino);
    print(pythia.particleData.findParticle(heavy_neutrino)->sizeChannels());
    print(tau);
    print(pythia.particleData.mWidth(heavy_neutrino));


    HepMC::IO_GenEvent hepmc_file(file_name, std::ios::in);
    int total = 0;
    int good = 0;
    auto analysis = mapp::analysis();

    while (auto* hepmc_event = hepmc_file.read_next_event()) {
        ++total;
        pythia.event.reset();
        pythia.event.append(retrive_neutrino(hepmc_event, tau));
        if (!pythia.next()) {
            print("Pythia encountered a problem");
            continue;
        }
        if (debug) pythia.event.list(true);


        for (auto line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            auto vertex = particle.vProd();
            if (particle.chargeType() == 0) continue;
            if (!analysis.is_inside(get_point(particle))) continue;
            print("Hooray!", vertex.pAbs());
            ++good;
            break;
        }
        delete hepmc_event;
        hepmc_file >> hepmc_event;
    }
    auto fraction = double(good) / total;
    print("Fraction", fraction);
    print("sigma", sigma * fraction, "mb");
    return 0;
}

int main()
{
    return read_hepmc("neutrino_0.500000.hep");
    return write_branching_fractions();
    return write_hepmc(1.);
    return write_hepmcs();
}

