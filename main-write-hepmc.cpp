#include "pythia.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    hnl::write_hepmc(arguments.empty() ? .5 : hnl::convert(arguments.front()));
}
