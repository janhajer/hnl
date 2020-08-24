#include "string.hh"
#include "read-file.hh"
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

using namespace hnl;

std::vector<std::string> split_line(std::string const& line) {
    std::vector<std::string> strings;
    boost::split(strings, line, [](char c) {
        return c == ' ';
    }, boost::token_compress_on);
    return strings;
}

using Res = boost::optional<double>;
Res do_find(std::vector<std::string> const& vector, int which) {
    switch (which) {
        case 0 : return vector.size() == 6 && vector[0] == "HNLs" && vector[1] == "with" && vector[2] == "m" && vector[3] == "=" && vector[5] == "GeV" ? Res(to_double(vector[4])) : boost::none;
        case 1 : return vector.size() > 14 && vector[0] == "Events" && vector[1] == "were" && vector[2] == "generated" && vector[3] == "with" && vector[4] == "U^2" ? Res(to_double(vector[13])) : boost::none;
        case 2 : return vector.size() >= 3 && vector[0] == "In" && vector[1] == "MAPP" && vector[3] == "mb" ? Res(to_double(vector[2])) : boost::none;
        default : return boost::none;
    }
}

auto find_in_file(std::vector<std::string> const& lines) {
    int missing = 0;
    std::map<double, std::map<double, double>> result;
    std::array<double, 3> array;
    for (auto const& line : lines) {
        auto res = do_find(split_line(boost::trim_copy_if(line, boost::is_any_of("\t "))), missing);
        if(res) {
            array[missing] = *res;
        ++missing;
        }
        if (missing == 2) {
            print(array[0], array[1], array[2]);
            result[array[0]][array[1]] = array[2];
            missing = 0;
        }
    }
    return result;
}

int main(int argc, char** argv) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    save_result(find_in_file(import_file(arguments.empty() ? "test.out" : arguments.front())), arguments.empty() ? "test.out" : arguments.front());
}
