#include "read-lhe.hh"

int main(int argc, char** argv) {
    using namespace hnl;
    std::vector<std::string> arguments(argv + 1, argv + argc);
    print(read_lhe(arguments.empty() ? "test.lhe" : arguments.front(), 1));
}
