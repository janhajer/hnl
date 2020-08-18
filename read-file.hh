#pragma once

#include <bits/stdc++.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "string.hh"
#include "id.hh"

namespace hnl {

namespace{
    bool const debug_2 = false;
}

auto split_line(std::string const& line) {
    std::vector<std::string> strings;
    boost::split(strings, line, [](char c) {
        return c == ' ';
    }, boost::token_compress_on);
    return strings;
}

template<typename Predicate>
auto find_in_file_copy(std::vector<std::string>& lines, int pos, Predicate predicate) {
    auto found = boost::range::find_if(lines, [&predicate](auto & line) {
        boost::trim_if(line, boost::is_any_of("\t "));
        return predicate(split_line(line));
    });
    return found == lines.end() ? "value not found"s : split_line(*found).at(pos);
}

// template<typename Predicate>
// auto read_file(boost::filesystem::path const& path, int pos, Predicate predicate) {
//     std::ifstream file(path.string());
//     std::vector<std::string> lines;
//     std::copy(std::istream_iterator<Line>(file), std::istream_iterator<Line>(), std::back_inserter(lines));
//     auto found = boost::range::find_if(lines, [&predicate](auto & line) {
//         boost::trim_if(line, boost::is_any_of("\t "));
//         return predicate(split_line(line));
//     });
//     return found == lines.end() ? "value not found"s : split_line(*found).at(pos);
// }

auto find_sigma(std::vector<std::string>& lines) {
    if (debug_2) print("find sigma");
    return find_in_file_copy(lines, 1, [](auto const & strings)  {
        return /*strings.size() == 3 &&*/ strings.at(0) == "sigma" && strings.at(2) == "mb";
    });
}

auto find_mass(std::vector<std::string>& lines) {
    if (debug_2) print("find mass");
    return find_in_file_copy(lines, 1, [](auto const & strings)  {
        return strings.size() == 3 && strings.at(0) == "mass" && strings.at(2) == "GeV";
    });
}

auto find_coupling(std::vector<std::string>& lines, int heavy, int light) {
    if (debug_2) print("find coupling");
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
    while (std::fgets(string, sizeof(string), file)) lines.emplace_back(string);
    for (auto & line : lines) line.pop_back();
    return lines;
}

struct Line {
    friend std::istream& operator>> (std::istream& stream, Line& line) {
        std::getline(stream, line.string);
        return stream;
    }
    operator std::string() const {
        return string;
    }
private:
    std::string string;
};

auto import_head(boost::filesystem::path const& path, int number) {
    std::ifstream file(path.string());
    std::vector<std::string> lines;
    std::copy_n(std::istream_iterator<Line> (file), number, std::back_inserter(lines));
    return lines;
}
auto import_tail(boost::filesystem::path const& path, int number) {
    FILE* fp = std::fopen(path.string().c_str(), "r");
    auto back = tail(fp, number);
    fclose(fp); ;
    return back;
}


struct Meta {
    double mass = 0;
    double sigma = 0;
    std::map<int, std::map<int, double>> couplings;
};

boost::optional<Meta> meta_info(boost::filesystem::path const& path) {
    auto lines = import_head(path, 100) + import_tail(path, 100);
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

auto find_mass_lhe(std::vector<std::string>& lines) {
    return find_in_file_copy(lines, 1, [](auto const & strings) {
        return strings.size() > 2 && strings.at(0) == std::to_string(heavy_neutrino) && strings.at(2) == "#" && strings.at(3) == "mn1";
    });
}

auto find_sigma_lhe(std::vector<std::string>& lines) {
    return find_in_file_copy(lines, 5, [](auto const & strings) {
        return strings.size() > 4 && strings.at(0) == "#" && strings.at(1) == "Integrated" && strings.at(2) == "weight" && strings.at(3) == "(pb)" && strings.at(4) == ":";
    });
}

std::string get_param_heavy(int heavy){
    switch (heavy){
        case 9900012 : return "n1";
        case 9900014 : return "n2";
        case 9900016 : return "n3";
        default : print("not a heavy neutrino");
    }
    return "";
}

std::string get_param_light(int light){
    switch (light){
        case 12 : return "e";
        case 14 : return "mu";
        case 16 : return "ta";
        default : print("not a heavy neutrino");
    }
    return "";
}

std::string get_param(int heavy, int light){
    return "v" + get_param_light(light) + get_param_heavy(heavy);
}

std::string find_coupling_lhe(std::vector<std::string>& lines, int heavy, int light, int pos) {
    std::string name = get_param(heavy, light);
    return find_in_file_copy(lines, 1, [&name, pos](auto const & strings) noexcept {
        return strings.size() > 3 && strings.at(0) == std::to_string(pos) && strings.at(2) == "#" && strings.at(3) == name;
    });
}

boost::optional<Meta> meta_info_lhe(boost::filesystem::path const& path) {
    auto lines = import_head(path, 500);
    Meta meta;
    meta.mass = to_double(find_mass_lhe(lines));
    if (meta.mass <= 0) return boost::none;
    meta.sigma = to_double(find_sigma_lhe(lines));
    if (meta.sigma <= 0) return boost::none;
    int pos = 0;
    for (auto light : light_neutrinos()) for (auto heavy : heavy_neutral_leptons()) meta.couplings[heavy][light] = to_double(find_coupling_lhe(lines, heavy, light, ++pos));
    if (meta.couplings.empty()) return boost::none;
    print("Meta info",meta.mass, meta.sigma);
    return meta;
}

}
