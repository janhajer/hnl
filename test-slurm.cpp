#include "read-hepmc.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    hnl::hepmc::extract_metas(arguments.empty() ? "." : arguments.front());
}
