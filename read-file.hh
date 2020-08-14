#pragma once

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "id.hh"
#include "string.hh"

namespace hnl {

// namespace {
//
// constexpr bool debug = false;
//
// }

struct Line {
    friend std::istream& operator>>(std::istream& stream, Line& line) noexcept {
        std::getline(stream, line.string);
        return stream;
    }
    operator std::string() const noexcept {
        return string;
    }
private:
    std::string string;
};

auto split_line(std::string const& line) noexcept {
    std::vector<std::string> strings;
    boost::split(strings, line, [](char c) noexcept {
        return c == ' ';
    }, boost::token_compress_on);
    return strings;
}

auto import_file(boost::filesystem::path const& path) noexcept {
    std::ifstream file(path.string());
    std::vector<std::string> lines;
    std::copy(std::istream_iterator<Line>(file), std::istream_iterator<Line>(), std::back_inserter(lines));
    return lines;
}

template<typename Predicate>
auto read_file(std::vector<std::string>& lines, int pos, Predicate predicate) noexcept {
    auto found = boost::range::find_if(lines, [&predicate](auto & line) noexcept {
        boost::trim_if(line, boost::is_any_of("\t "));
        return predicate(split_line(line));
    });
    return found == lines.end() ? "value not found"s : split_line(*found).at(pos);
}

template<typename Predicate>
auto read_file(boost::filesystem::path const& path, int pos, Predicate predicate) noexcept {
    std::ifstream file(path.string());
    std::vector<std::string> lines;
    std::copy(std::istream_iterator<Line>(file), std::istream_iterator<Line>(), std::back_inserter(lines));
    auto found = boost::range::find_if(lines, [&predicate](auto & line) noexcept {
        boost::trim_if(line, boost::is_any_of("\t "));
        return predicate(split_line(line));
    });
    return found == lines.end() ? "value not found"s : split_line(*found).at(pos);
}

auto find_sigma(std::vector<std::string>& path) {
    if (debug) print("find sigma");
    return read_file(path, 1, [](auto const & strings)  {
        return strings.size() == 3 && strings.at(0) == "sigma" && strings.at(2) == "mb";
    });
}

auto find_mass(std::vector<std::string>& path) {
    if (debug) print("find mass");
    return read_file(path, 1, [](auto const & strings)  {
        return strings.size() == 3 && strings.at(0) == "mass" && strings.at(2) == "GeV";
    });
}

auto find_coupling(std::vector<std::string>& path, int heavy, int light) {
    if (debug) print("find coupling");
    return read_file(path, 3, [&](auto const & strings)  {
        return strings.size() == 4 && strings.at(0) == "coupling" && strings.at(1) == std::to_string(heavy) && strings.at(2) == std::to_string(light);
    });
}



}
