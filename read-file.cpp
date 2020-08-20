#include <bits/stdc++.h>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <boost/filesystem.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>

#include "read-file.hh"
#include "string.hh"

namespace hnl {

namespace {
bool const debug = false;
}

std::vector<std::string> split_line(std::string const& line) {
    std::vector<std::string> strings;
    boost::split(strings, line, [](char c) {
        return c == ' ';
    }, boost::token_compress_on);
    return strings;
}

std::string find_in_file_copy(std::vector<std::string>& lines, int pos, std::function<bool(std::vector<std::string>)> const& predicate) {
    auto found = boost::range::find_if(lines, [&predicate](auto & line) {
        boost::trim_if(line, boost::is_any_of("\t "));
        return predicate(split_line(line));
    });
    return found == lines.end() ? "value not found" : split_line(*found).at(pos);
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
    for (auto& line : lines) line.pop_back();
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

std::vector<std::string> import_head(boost::filesystem::path const& path, int number) {
    std::ifstream file(path.string(), std::ios_base::in | std::ios_base::binary);

    std::vector<std::string> lines;
    if (path.extension().string() == ".gz") {
        boost::iostreams::filtering_streambuf<boost::iostreams::input> buffer;
        buffer.push(boost::iostreams::gzip_decompressor());
        buffer.push(file);
        std::istream instream(&buffer);
        std::copy_n(std::istream_iterator<Line>(instream), number, std::back_inserter(lines));
    } else std::copy_n(std::istream_iterator<Line>(file), number, std::back_inserter(lines));
    return lines;
}

std::vector<std::string> import_tail(boost::filesystem::path const& path, int number) {
    // TODO move from C to C++
    FILE* fp = std::fopen(path.string().c_str(), "r");
    auto back = tail(fp, number);
    fclose(fp); ;
    return back;
}

}
