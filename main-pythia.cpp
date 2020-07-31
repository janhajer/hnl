#include "generic.hh"
#include "ResonanceWidths.hh"

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
// #include "HepMC/Units.h"

namespace
{

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

namespace neutrino
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
    using namespace neutrino;
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
    return id == 0 ? 1 : 0;
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
//     pythia.particleData.m0(1, 0.0048);
//     pythia.particleData.m0(2, 0.0023);
//     pythia.particleData.m0(3, 0.095);
//     pythia.particleData.m0(4, 1.275);
//     pythia.particleData.m0(5, 4.180);
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

void get_neutrino(double mass)
{
    using namespace neutrino;

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

}

void get_neutrinos()
{
    using namespace neutrino;
    Loop loop(.5, 10);
    for (auto step = 0; step <= loop.steps; ++step) {
        auto mass = loop.mass(6., step);
        std::ofstream ofstream("neutrino_" + std::to_string(mass) + ".txt");
        std::streambuf* streambuf = std::cout.rdbuf();
        std::cout.rdbuf(ofstream.rdbuf());
        get_neutrino(mass);
        std::cout.rdbuf(streambuf);
    }
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
//     pythia.readString("ResonanceWidths:minWidth = 1E-15");
}

auto has_neutrino = [](auto const& channel)
{
    return channel.product(0) == heavy_neutrino || channel.product(1) == heavy_neutrino || channel.product(2) == heavy_neutrino || channel.product(3) == heavy_neutrino || channel.product(4) == heavy_neutrino;
};

template<typename Data>
void save_data(Data& result, neutrino::Loop const& loop, double mass, int source)
{
    using namespace neutrino;
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

void get_branching_fractions()
{
    using Result = std::map<int, std::map<std::tuple<int, int, int, int, int>, std::map<int, double>>>;
    using namespace neutrino;

//     std::vector<int> sources{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541, 443, 553};
//     std::vector<int> sources{431, 411, 421};
//     std::vector<int> sources{511, 521, 531, 541};
    std::vector<int> sources{511, 521, 531};
//     std::vector<int> sources{443, 553};
//     std::vector<int> sources{heavy_neutrino};
    Result result;
    for (auto source : sources) {
        Loop loop(.1, 100);
        double mass_max = 5;
        for (auto step = 0; step <= loop.steps; ++step) {
            print(step, "of", loop.steps);
            Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
            set_pythia(pythia);
            if (!is_neutrino(source)) mass_max = pythia.particleData.m0(source) * 1.01;
            pythia.particleData.m0(heavy_neutrino, loop.mass(mass_max, step));
            is_neutrino(source) ? pythia.setResonancePtr(new NeutrinoResonance(pythia, neutrino_coupling_2, source)) : pythia.setResonancePtr(new MesonResonance(pythia, neutrino_coupling_2, source));
            pythia.particleData.findParticle(source)->rescaleBR();
            pythia.init();
            auto const& particle = *pythia.particleData.findParticle(source);
            for (auto pos = 0; pos < particle.sizeChannels(); ++pos) {
                auto channel = particle.channel(pos);
                auto ratio = channel.bRatio();
                if (ratio > 0. && (has_neutrino(channel) || is_neutrino(source)))
                    result[source][ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}][step] = ratio;
            }
        }
        save_data(result, loop, mass_max, source);
    }
}

void read_neutrinos_1()
{
    std::string file = "weakbosons";
    Pythia8::Pythia pythia;
    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF =" + file + ".lhe");
    pythia.init();

    for (int event = 0; ; ++event) {
        if (!pythia.next()) {
            if (pythia.info.atEndOfFile()) break;
            continue;
        }
        for (int line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            neutrino::print("next line");
            if (std::abs(particle.id()) != heavy_neutrino) continue;
            neutrino::print("pt", particle.pT(), "eta", particle.eta(), "phi", particle.phi());
        }
    }
    pythia.stat();
}

using namespace neutrino;
auto is_good_event = [](HepMC::GenEvent* const event) -> bool
{
    for (auto iterator = event->particles_begin(); iterator != event->particles_end(); ++iterator) {
        auto* particle = *iterator;
        if (particle->pdg_id() == heavy_neutrino) {
            auto vProdSave = particle->production_vertex()->position();
            auto pSave = particle->momentum();
            auto tauSave = 1.;
            auto mSave = pSave.m();

            double xDec = vProdSave.x() + tauSave * pSave.px() / mSave;
            double yDec = vProdSave.y() + tauSave * pSave.py() / mSave;
            double zDec = vProdSave.z() + tauSave * pSave.pz() / mSave;
            double tDec = vProdSave.t() + tauSave * pSave.e() / mSave;
            auto dec = HepMC::FourVector(xDec, yDec, zDec, tDec);


            return true;
        }
    }
    return false;
};

void read_neutrinos_2(double)
{
    HepMC::IO_GenEvent ascii_in("test", std::ios::in);
    int total = 0;
    int good = 0;
    auto event = ascii_in.read_next_event();
    while (event) {
        total++;
        if (is_good_event(event)) {
            ++good;
        }
        delete event;
        ascii_in >> event;
    }
    double fraction = double(good) / total;
    print("Fraction", fraction);
}

int main(int, char* [])
{
//     get_neutrino(1.);
//     read_neutrinos_1();
//     read_neutrinos_2();
    get_neutrinos();
//     get_branching_fractions();
}

