#include <iostream>
#include <fstream>
#include <vector>
#include <string>

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    std::ofstream file(arguments.empty() ? "testing-slurm.dat" : arguments.at(0) + ".dat");
    file << "Writing this to a file.\n";
    file << (arguments.empty() ? "testing-slurm.txt" : arguments.at(0) + ".txt");
}
