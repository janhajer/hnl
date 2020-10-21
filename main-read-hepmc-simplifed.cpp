#include "read-hepmc.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    hnl::hepmc::read_simplified(
        arguments.size() < 1 ? "0.500000.hep" : arguments.at(0),
        arguments.size() < 2 ? 1 : hnl::to_double(arguments.at(1))
    );
}
