#include "write-hepmc.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    hnl::calculate_sigma(arguments.empty() ? 2. : hnl::to_double(arguments.front()));
}
