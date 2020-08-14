#pragma once

#include <bits/stdc++.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "string.hh"

namespace hnl {

auto split_line(std::string const& line) noexcept {
    std::vector<std::string> strings;
    boost::split(strings, line, [](char c) noexcept {
        return c == ' ';
    }, boost::token_compress_on);
    return strings;
}

template<typename Predicate>
auto find_in_file_copy(std::vector<std::string>& lines, int pos, Predicate predicate) noexcept {
    auto found = boost::range::find_if(lines, [&predicate](auto & line) noexcept {
        boost::trim_if(line, boost::is_any_of("\t "));
        return predicate(split_line(line));
    });
    return found == lines.end() ? "value not found"s : split_line(*found).at(pos);
}

// template<typename Predicate>
// auto read_file(boost::filesystem::path const& path, int pos, Predicate predicate) noexcept {
//     std::ifstream file(path.string());
//     std::vector<std::string> lines;
//     std::copy(std::istream_iterator<Line>(file), std::istream_iterator<Line>(), std::back_inserter(lines));
//     auto found = boost::range::find_if(lines, [&predicate](auto & line) noexcept {
//         boost::trim_if(line, boost::is_any_of("\t "));
//         return predicate(split_line(line));
//     });
//     return found == lines.end() ? "value not found"s : split_line(*found).at(pos);
// }

auto find_sigma(std::vector<std::string>& lines) {
    if (debug) print("find sigma");
    return find_in_file_copy(lines, 1, [](auto const & strings)  {
        return strings.size() == 3 && strings.at(0) == "sigma" && strings.at(2) == "mb";
    });
}

auto find_mass(std::vector<std::string>& lines) {
    if (debug) print("find mass");
    return find_in_file_copy(lines, 1, [](auto const & strings)  {
        return strings.size() == 3 && strings.at(0) == "mass" && strings.at(2) == "GeV";
    });
}

auto find_coupling(std::vector<std::string>& lines, int heavy, int light) {
    if (debug) print("find coupling");
    return find_in_file_copy(lines, 3, [&](auto const & strings)  {
        return strings.size() == 4 && strings.at(0) == "coupling" && strings.at(1) == std::to_string(heavy) && strings.at(2) == std::to_string(light);
    });
}


std::vector<std::string> tail(FILE* file, int n) {
    std::vector<std::string> lines;
    int count = 0;
    char string[2 * 100];
    if (std::fseek(file, 0, SEEK_END)) return lines;
    auto pos = std::ftell(file);
    while (pos) {
        if (std::fseek(file, --pos, SEEK_SET)) return lines;
        if (std::fgetc(file) == '\n') if (count++ == n) break;
    }
    while (std::fgets(string, sizeof(string), file)) {
        std::string  sstring;
        std::getline(string,sstring);
        lines.emplace_back(sstring);
    }
    return lines;
}

struct Line {
    friend std::istream& operator>> (std::istream& stream, Line& line) noexcept {
        std::getline(stream, line.string);
        return stream;
    }
    operator std::string() const noexcept {
        return string;
    }
private:
    std::string string;
};

auto import_lines(boost::filesystem::path const& path) noexcept {
    std::ifstream file(path.string());
    std::vector<std::string> lines;
    std::copy_n(std::istream_iterator<Line> (file), 100, std::back_inserter(lines));
    print(lines);
    FILE* fp = std::fopen(path.string().c_str(), "r");
    auto back = tail(fp, 100);
    fclose(fp); ;
    print(back);
    return lines + back;
}


struct Meta {
    double mass = 0;
    double sigma = 0;
    std::map<int, std::map<int, double>> couplings;
};

boost::optional<Meta> meta_info(boost::filesystem::path const& path) {
    auto lines = import_lines(path);
    Meta meta;
    meta.mass = to_double(find_mass(lines));
    if (meta.mass <= 0) return boost::none;
    meta.sigma = to_double(find_sigma(lines));
    if (meta.sigma <= 0) return boost::none;
    for (auto heavy : heavy_neutral_leptons()) for (auto light : light_neutrinos()) meta.couplings[heavy][light] = to_double(find_coupling(lines, heavy, light));
    if (meta.couplings.empty()) return boost::none;
    print("Meta info",meta.mass, meta.sigma);
    return meta;
}

}
