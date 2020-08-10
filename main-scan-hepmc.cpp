#include "pythia.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    hnl::scan_hepmc(arguments.empty() ? "0.500000.hep" : arguments.front());
}
