#include "Pythia8/Pythia.h"
#include "units.hh"
#include "id.hh"
#include "Decayer.hh"
#include "ResonanceWidths.hh"
#include "ResonanceWidth.hh"
// #include <boost/units/pow.hpp>
#include <iostream>
#include <map>
#include <tuple>

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
    UserHook(double eta_min_, double eta_max_) : eta_min(eta_min_), eta_max(eta_max_)  {}
    virtual bool canVetoStep() override
    {
        return true;
    }
    virtual bool doVetoResonanceDecays(Pythia8::Event& event) override
    {
        for (auto line = 0; line < event.size(); ++line) {
            auto particle = event[line];
            if (std::abs(particle.id()) > 200) print(line, particle.id());
            if (std::abs(particle.id()) == to_underlying(Id::neutrino2)) {
//                 && std::abs(particle.eta()) > eta_min && std::abs(particle.eta()) < eta_max)
                print("Found !!");
                return false;
            }
        }
        print("Not found");
        return true;
    }
private:
    double eta_min;
    double eta_max;
};

}

auto get_hundredth2(int id)
{
    return id > 999 ? 100 * ((id / 100) % 10) + 10 * ((id / 10) % 10) + id % 10 : id;
}

struct BRatio {
    int onMode;
    double bRatio;
    int meMode;
};

int main(int argc, char* argv[])
{
    using namespace neutrino;

    Pythia8::Pythia pythia;
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 14000.");

    pythia.readString("Bottomonium:all = on");
// pythia.readString("Charmonium:all = on");

    std::map<std::tuple<int, int, int,int , int>, BRatio> tmp;
    auto const* part = pythia.particleData.findParticle(521);
    for (auto n = 0; n < part->sizeChannels(); ++n) {
        auto channel = part->channel(n);
        print(channel.product(0), channel.product(1), channel.product(2), channel.onMode(), channel.bRatio(), channel.meMode());
        tmp[ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}] = {channel.onMode(), channel.bRatio(), channel.meMode()};
    }


    std::map<std::tuple<int, int, int, int, int>, std::map<int, double>> result;
//     std::vector<std::vector<double>> table;
    int iterations = 10;
    int min = 0;
    double max = 5;
    double offset = 0.;
    for (auto i = min; i <= iterations; ++i) {
        double mass = max * i / iterations + offset;
        pythia.particleData.m0(9900014, mass);
        auto* ptr = new MesonResonance(pythia.settings, &pythia.rndm, 1, 521);
        ptr->init(&pythia.info, &pythia.settings, &pythia.particleData, &pythia.couplings);
        pythia.particleData.initWidths({ptr});
        pythia.init();
//         std::map<double> row{mass};
        std::map<double, double> row;
        auto * part = pythia.particleData.findParticle(521);
        for(auto const& row : tmp) part->addChannel(row.second.onMode, row.second.bRatio, row.second.meMode, std::get<0>(row.first), std::get<1>(row.first), std::get<2>(row.first), std::get<3>(row.first), std::get<4>(row.first));
        for (auto n = 0; n < part->sizeChannels(); ++n) {
            auto channel = part->channel(n);
            auto ratio = channel.bRatio();
            print(channel.product(0), channel.product(1), channel.product(2), ratio);
            if (ratio > 0.) result[ {channel.product(0), channel.product(1), channel.product(2), channel.product(3), channel.product(4)}][i] = ratio;
        }
//         table.emplace_back(row);
    }
//     print("res", table);
    std::ofstream output_file("test.dat");
    output_file << 0 << '\t' << 1 << '\t' << 2 << '\t' << 3 << '\t' << 4;
    for (auto i = 0; i <= iterations; ++i) output_file << std::scientific << '\t' << max* i / iterations + offset;
    output_file << '\n';

    for (auto& row : result) {
        output_file << std::scientific << std::get<0>(row.first) << '\t' << std::get<1>(row.first) << '\t' << std::get<2>(row.first) << '\t' << std::get<3>(row.first) << '\t' << std::get<4>(row.first) << '\t';
        for (auto i = 0; i < iterations; ++i) {
            output_file << std::scientific << row.second[i] << '\t';
//         std::copy(row.cbegin(), row.cend(), std::ostream_iterator<double>(output_file << std::scientific, '\t'));
        }
        output_file << '\n';

    }

}

int main_test(int argc, char* argv[])
{
    using namespace neutrino;

    Pythia8::Pythia pythia;
    pythia.particleData.m0(9900014, 1);
    pythia.readString("9900014:addChannel = 1 0. 0 13 -13 12");
    pythia.readString("9900014:addChannel = 1 0. 0 13 -13 12");
    pythia.readString("9900014:addChannel = 1 0. 0 13 -13 14");
    pythia.readString("9900014:addChannel = 1 0. 0 13 -13 16");

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

    pythia.particleData.initWidths({new NeutrinoResonanceWidth(9900014)});

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
