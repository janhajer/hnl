#include "container.hh"
#include "read-file.hh"
#include "read-lhe.hh"

int main(int argc, char** argv) {
    using namespace hnl;
    std::vector<std::string> arguments(argv + 1, argv + argc);
    auto meta = meta_info_lhe(arguments.empty() ? "test.lhe" : arguments.front());
    if(meta) print(meta->mass, meta->sigma,meta->couplings);
}
