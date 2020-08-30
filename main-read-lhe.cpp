#include "read-lhe.hh"

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    hnl::print(hnl::lhe::read(arguments.empty() ? "test.lhe" : arguments.front(), 1));
}
