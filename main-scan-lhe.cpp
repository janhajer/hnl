#include "read-lhe.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    hnl::scan_lhe(arguments.empty() ? "test.lhe.gz" : arguments.front());
}
