#include "pythia.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    hnl::write_sigma_hepmc(arguments.empty() ? 2. : hnl::to_double(arguments.front()));
}
