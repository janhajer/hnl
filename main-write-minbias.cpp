#include "write-hepmc.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    hnl::write_minimum_bias(arguments.empty() ? .5 : hnl::to_double(arguments.front()));
}
