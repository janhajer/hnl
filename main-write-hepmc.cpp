#include "write-hepmc.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    double mass = 0., coupling = 0.;
    int events = 1000, seed = -1;
    if (arguments.size() < 3 || arguments.size() > 5
        || !hnl::parse(arguments.at(0), mass) || mass <= 0.
        || !hnl::parse(arguments.at(1), coupling) || coupling <= 0.
        || hnl::hepmc::parents(arguments.at(2)).empty()
        || (arguments.size() > 3 && (!hnl::parse(arguments.at(3), events) || events <= 0))
        || (arguments.size() > 4 && (!hnl::parse(arguments.at(4), seed) || seed < 0 || seed > 900000000))) {
        hnl::print("usage: write-hepmc <m_N [GeV] > 0> <|U_mu|^2 > 0> <parent flavour: b, c, s or pi> [events > 0, default 1000] [seed]");
        return 1;
    }
    return hnl::hepmc::write_hepmc(mass, coupling, arguments.at(2), events, seed) ? 0 : 1;
}
