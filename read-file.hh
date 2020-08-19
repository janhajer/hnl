#pragma once

#include <bits/stdc++.h>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm/find_if.hpp>

#include "string.hh"

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
    std::ifstream file(path.string(), std::ios_base::in | std::ios_base::binary);

    boost::iostreams::filtering_streambuf<boost::iostreams::input> buffer;
    buffer.push(boost::iostreams::gzip_decompressor());
    buffer.push(file);
    std::istream instream(&buffer);


    std::vector<std::string> lines;
//     std::copy_n(std::istream_iterator<Line>(file), number, std::back_inserter(lines));
    std::copy_n(std::istream_iterator<Line>(instream), number, std::back_inserter(lines));
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

}
