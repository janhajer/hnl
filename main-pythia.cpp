#include "Pythia8/Pythia.h"
#include "units.hh"
#include "id.hh"
#include "Decayer.hh"
#include "ResonanceWidths.hh"
#include "ResonanceWidth.hh"
#include "Pythia8Plugins/HepMC2.h"
#include "HepMC/Units.h"
// #include <boost/units/pow.hpp>
#include <iostream>
#include <map>
#include <tuple>



namespace
{

const int heavy_neutrino = 9900012;

}


namespace neutrino
{

std::vector<std::string> three_body_decay(int in, std::vector<std::pair<int, double>> out)
{
    return transform(out, [in](auto const & pair) {
        return std::to_string(in) + ":addChannel = 1 " + std::to_string(pair.second) + " 0 -13 9900014 " + (pair.first == 0 ? "" : std::to_string(pair.first));
    }) + transform(out, [in](auto const & pair) {
        return std::to_string(in) + ":addChannel = 1 " + "0." + " 0 -13 14 " + (pair.first == 0 ? "" : std::to_string(pair.first));
    });
}

std::vector<std::string> four_body_decay(int in, std::vector<std::tuple<int, int, double>> out)
{
    return transform(out, [in](auto const & tuple) {
        return std::to_string(in) + ":addChannel = 1 " + std::to_string(std::get<2>(tuple)) + " 0 -13 9900014 " + std::to_string(std::get<0>(tuple)) + " " + std::to_string(std::get<1>(tuple));
    }) + transform(out, [in](auto const & tuple) {
        return std::to_string(in) + ":addChannel = 1 " + "0." + " 0 -13 14 " + std::to_string(std::get<0>(tuple)) + " " + std::to_string(std::get<1>(tuple));
    });
}

struct UserHook : public Pythia8::UserHooks {
    UserHook(double eta_min_ = 0, double eta_max_ = 0) : eta_min(eta_min_), eta_max(eta_max_)  {}
    virtual bool canVetoPartonLevel() override
    {
        return true;
    }
    virtual int numberVetoStep() override
    {
        return 100;
    }
    virtual bool doVetoPartonLevel(Pythia8::Event const& event) override
    {
        print("doVetoStep", event.size());
        for (auto line = 0; line < event.size(); ++line) {
            auto const& particle = event[line];
            auto absid = std::abs(particle.id());
            if (absid == heavy_neutrino) {
                print("Success in", "line", line, "from mother", event[particle.mother1()].id(), "with decay vertex", particle.vDec().pAbs());
                return false;
            }
        }
        print("Not found");
        return false;
    }
private:
    double eta_min;
    double eta_max;
};


auto log_scale(double min, double max, int bin, int bin_number) noexcept
{
    auto log_min = std::log10(min);
    auto log_max = std::log10(max);
    auto bin_width = (log_max - log_min) / bin_number;
//     return transform(irange(bin_number), [&](auto bin) noexcept {
    return std::pow(10, log_min + bin * bin_width);
//     });
}


auto lin_scale(double min, double max, int bin, int bin_number) noexcept
{
    auto bin_width = (max - min) / bin_number;
    return min + bin * bin_width;
}

auto get_hundredth2(int id)
{
    return id > 999 ? 100 * ((id / 100) % 10) + 10 * ((id / 10) % 10) + id % 10 : id;
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
    return id_heavy == 9900012 && id_light == 12 ? 1 : 0;
};

auto get_resonances(Pythia8::Pythia& pythia, std::vector<int> const& mesons)
{
    std::vector<Pythia8::ResonanceWidths*> resonances;
    for (auto meson : mesons) resonances.emplace_back(new MesonResonance(pythia, neutrino_coupling_2, meson));
    return resonances;
}

struct Loop {
    double m_min = .1;
    int steps = 20;
    double mass(double max, int step)
    {
        return lin_scale(m_min, max, step, steps);
    }
};

}

void main_single(double mass)
{
    using namespace neutrino;

    HepMC::Pythia8ToHepMC pythia_to_hep;
    pythia_to_hep.set_store_proc();
    pythia_to_hep.set_store_xsec();
    pythia_to_hep.set_store_pdf();
    HepMC::IO_GenEvent io_event("neutrino_" + std::to_string(mass) + ".hep", std::ios::out);

    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000.");
//     pythia.readString("Bottomonium:all = on");
//     pythia.readString("Charmonium:all = on");
//     pythia.readString("Onia:all(3S1) = on");
    pythia.readString("SoftQCD:nonDiffractive = on");
//     pythia.readString("HardQCD:qqbar2ccbar  = on");
//     pythia.readString("HardQCD:gg2bbbar  = on");
//     pythia.readString("HardQCD:all = on");
//     pythia.readString("PhaseSpace:pTHatMin = .1");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Main:numberOfEvents = 1000");
    pythia.particleData.m0(heavy_neutrino, mass);

    std::vector<int> mesons{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541, 443, 553};
    auto resonances = transform(mesons, [&pythia](int meson) -> Pythia8::ResonanceWidths* {
        return new MesonResonance(pythia, neutrino_coupling_2, meson);
    });
//     auto resonances = get_resonances(pythia, mesons);
    pythia.particleData.initWidths(resonances);
    for (auto meson : mesons) pythia.particleData.findParticle(meson)->rescaleBR();

    pythia.init();

    int total = 0;
    int successfull = 0;
    while (successfull < pythia.mode("Main:numberOfEvents")) {
        if (!pythia.next()) continue;
        ++total;
        bool success = false;
        for (auto line = 0; line < pythia.event.size(); ++line) {
            auto const& particle = pythia.event[line];
            if (std::abs(particle.id()) != heavy_neutrino) continue;
            double eta = std::abs(particle.eta());
            if (!(1.4 < eta && eta < 3.5)) continue;
            print("Success event", successfull, "of", total, "in line", line, "from mother", pythia.event[particle.mother1()].id(), "with decay vertex", particle.vDec().pAbs(), "and", particle.eta(), particle.phi());
            success = true;
            break;

        }
        if (!success) continue;
        ++successfull;
        HepMC::GenEvent gen_event;
        pythia_to_hep.fill_next_event(pythia, &gen_event);
        io_event.write_event(&gen_event);
    }

    pythia.stat();

    print("\nwritten", successfull, "events of", total, "that is a fraction", double(successfull) / total);
    print("\ntherefore", pythia.info.sigmaGen(), pythia.info.sigmaGen() * successfull / total);

    io_event.write_comment("sigma " + std::to_string(pythia.info.sigmaGen() * successfull / total) + " mb");

}

void loop_single()
{
    using namespace neutrino;
    Loop loop;
    for (auto step = 0; step <= loop.steps; ++step) {
        auto mass = loop.mass(4., step);
        std::ofstream ofstream("neutrino_" + std::to_string(mass) + ".txt");
        std::streambuf* streambuf = std::cout.rdbuf();
        std::cout.rdbuf(ofstream.rdbuf());
        main_single(mass);
        std::cout.rdbuf(streambuf);
    }
}

void set_pythia(Pythia8::Pythia& pythia)
{
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000.");
    pythia.readString("Bottomonium:all = on");

    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Init:showProcesses  = off");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showMultipartonInteractions = off");
    pythia.readString("ResonanceWidths:minWidth = 1E-25");
}

auto has_neutrino = [](auto const& channel)
{
    return channel.product(0) == heavy_neutrino || channel.product(1) == heavy_neutrino || channel.product(2) == heavy_neutrino || channel.product(3) == heavy_neutrino || channel.product(4) == heavy_neutrino;
};


template<typename Data>
void save_data(Data& result, double mass, int meson)
{
    using namespace neutrino;
    Loop loop;
    std::ofstream output_file(std::to_string(meson) + ".dat");
    output_file << 0 << '\t' << 1 << '\t' << 2 << '\t' << 3 << '\t' << 4;
    for (auto step = 0; step <= loop.steps; ++step) output_file << std::scientific << '\t' << loop.mass(mass, step);
    output_file << '\n';
    for (auto& row : result[meson]) {
        output_file << std::scientific << std::get<0>(row.first) << '\t' << std::get<1>(row.first) << '\t' << std::get<2>(row.first) << '\t' << std::get<3>(row.first) << '\t' << std::get<4>(row.first) << '\t';
        for (auto step = 0; step <= loop.steps; ++step) output_file << std::scientific << row.second[step] << '\t';
        output_file << '\n';
    }
}

void main_loop()
{
    using namespace neutrino;

    std::vector<int> mesons{211, 130, 310, 321, 411, 421, 431, 511, 521, 531, 541, 443, 553};
//     std::vector<int> mesons{431,411,421};
//     std::vector<int> mesons{511, 521, 531, 541};
//     std::vector<int> mesons{443, 553};
    std::map<int, std::map<std::tuple<int, int, int, int, int>, std::map<int, double>>> result;
    for (auto meson : mesons) {
        Loop loop;
        double mass = 0;
        for (auto step = 0; step <= loop.steps; ++step) {
            print(step, "of", loop.steps);
            Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
            set_pythia(pythia);
            mass = pythia.particleData.m0(meson);
            pythia.particleData.m0(heavy_neutrino, loop.mass(mass, step));
            auto resonances = get_resonances(pythia, mesons);
            pythia.particleData.initWidths(resonances);
            for (auto* resonance : resonances) static_cast<MesonResonance*>(resonance)->AddMissingChannels(pythia.particleData);
            pythia.init();
//         for (auto meson : mesons) {
            auto const& particle = *pythia.particleData.findParticle(meson);
            for (auto pos = 0; pos < particle.sizeChannels(); ++pos) {
                auto channel = particle.channel(pos);
                auto ratio = channel.bRatio();
                if (ratio > 0. && has_neutrino(channel))
                    result[meson][ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}][step] = ratio;
            }
//         }
//     for (auto meson : mesons)
        }
        save_data(result, mass, meson);
    }
}

void main_read()
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

            double xDec = vProdSave.px() + tauSave * pSave.px() / mSave;
            double yDec = vProdSave.py() + tauSave * pSave.py() / mSave;
            double zDec = vProdSave.pz() + tauSave * pSave.pz() / mSave;
            double tDec = vProdSave.e()  + tauSave * pSave.e()  / mSave;
            auto dec = Vec4();


            return true;
        }
    }
    return false;
};

void main_read_2(double coupling)
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
//     main_single(1.);
//     main_read();
//     main_read_2();
//     loop_single();
    main_loop();
}


int main_test(int argc, char* argv[])
{
    using namespace neutrino;


    Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);
    pythia.particleData.m0(heavy_neutrino, 1);
    pythia.readString(std::to_string(heavy_neutrino) + ":addChannel = 1 0. 0 13 -13 12");
    pythia.readString(std::to_string(heavy_neutrino) + ":addChannel = 1 0. 0 13 -13 12");
    pythia.readString(std::to_string(heavy_neutrino) + ":addChannel = 1 0. 0 13 -13 14");
    pythia.readString(std::to_string(heavy_neutrino) + ":addChannel = 1 0. 0 13 -13 16");

    for (auto const& line : three_body_decay(411, {{0, 0.0004000}, {111, 0.0043000}, {113, 0.0028000}, {221, 0.0026000}, {223, 0.0028000}, {311, 0.0874000}, {-313, 0.0533000}, {-315, 0.0038000}, {331, 0.0005000}, {-10313, 0.0036000}})) pythia.readString(line);
    for (auto const& line : four_body_decay(411, {{311, 111, 0.0014000}, {-321, 211, 0.0027000}})) pythia.readString(line);

    for (auto const& line : three_body_decay(421, {{-211, 0.0034000}, {-213, 0.0022000}, {-321, 0.0340000}, {-323, 0.0214000}, {-325, 0.0015000}, {-10323, 0.0014000}})) pythia.readString(line);
    for (auto const& line : four_body_decay(421, {{311, -211, 0.0011000}, {-321, 111, 0.0006000}})) pythia.readString(line);

    for (auto const& line : three_body_decay(511, {{-211, 0.0001330}, { -213, 0.0002690}, { -411, 0.0207000}, { -413, 0.0570000}, { -415, 0.0023000}, { -10411, 0.0045000}, { -10413, 0.0052000}, { -20413, 0.0083000}})) pythia.readString(line);
    for (auto const& line : four_body_decay(511, {{-411, 111, 0.0010000}, {-413, 111, 0.0003000}, {-421, -211, 0.0020000}, {-423, -211, 0.0007000}, {1, -2, 0.0018920}, {-2, 1, 0.0018920}})) pythia.readString(line);

    for (auto const& line : three_body_decay(521, {{111, 0.0000720}, {113, 0.0001450}, {221, 0.0000840}, {223, 0.0001450}, {331, 0.0000840}, {-421, 0.0224000}, {-423, 0.0617000}, {-425, 0.0030000}, {-10421, 0.0049000}, {-10423, 0.0056000}, {-20423, 0.0090000}})) pythia.readString(line);
    for (auto const& line : four_body_decay(521, {{-411, 211, 0.0019000}, {-413, 211, 0.0006000}, {-421, 111, 0.0010000}, {-423, 111, 0.0003000}, {2, -2, 0.0019480}, {-2, 2, 0.0019480}})) pythia.readString(line);

    pythia.readFile(argc > 1 ? argv[1] : "mymain.cmnd");

    neutrino::Decayer decayer(&pythia.particleData, &pythia.rndm);
    std::vector<int> particles{411, 421, 511, 521};
    pythia.setDecayPtr(&decayer, particles);

    pythia.particleData.initWidths({new NeutrinoResonanceWidth(heavy_neutrino)});

//     UserHook user_hook(1.4, 3.5);
//     pythia.setUserHooksPtr(&user_hook);

    pythia.init();

    int total = 0;
    int successfull = 0;
    while (successfull < pythia.mode("Main:numberOfEvents")) {
        ++total;
        print("start");
        if (!pythia.next()) continue;
        print("keep going");
        ++successfull;
        for (auto line = 0; line < pythia.event.size(); ++line) {
            auto particle = pythia.event[line];
            auto absid = std::abs(particle.id());
            if (absid == to_underlying(Id::neutrino2)) print("Success in even ", total, " of ", successfull, " line ", line, " mother ", pythia.event[particle.mother1()].id(), " pos ", particle.vDec().pAbs());
        }
    }

    pythia.stat();

    return 0;
}
