#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

#include "container.hh"
#include "string.hh"
#include "read-file.hh"

using namespace hnl;

namespace {

const bool debug = false;

}

using Res = boost::optional<double>;

Res extract_from_line(std::vector<std::string> const& vector, int which) {
    switch (which) {
        case 0 : return vector.size() == 6 && vector[0] == "HNLs" && vector[1] == "with" && vector[2] == "m" && vector[3] == "=" && vector[5] == "GeV" ? Res(to_double(vector[4])) : boost::none;
        case 1 : return vector.size() > 14 && vector[0] == "Events" && vector[1] == "were" && vector[2] == "generated" && vector[3] == "with" && vector[4] == "U^2" ? Res(to_double(vector[13])) : boost::none;
        case 2 : return vector.size() >= 3 && vector[0] == "In" && vector[1] == "MAPP" && vector[3] == "mb" ? Res(to_double(vector[2])) : boost::none;
        default : return boost::none;
    }
}

auto find_in_file(std::vector<std::string> const& lines) {
    int missing = 0;
    Result result;
    std::array<double, 3> array;
    for (auto const& line : lines) {
        auto found = extract_from_line(split_line(boost::trim_copy_if(line, boost::is_any_of("\t "))), missing);
        if (!found) continue;
        array[missing++] = *found;
        if (missing == 3) {
            if (debug) print(array[0], array[1], array[2]);
            result[array[0]][array[1]] = array[2];
            missing = 0;
        }
    }
    return result;
}

auto read_out(boost::filesystem::path const& name) {
    return find_in_file(import_file(name));
}

int main_2(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    auto name = arguments.empty() ? "test.out" : arguments.front();
    save(read_out(name), name);
    return 0;
}

auto reads_out(boost::filesystem::path const& path_name) {
    Result result;
    for (auto const& file : files(path_name)) if (file.path().extension().string() == ".out") result += read_out(file.path());
    save(result);
}

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    reads_out(arguments.empty() ? "." : arguments.front());
}
