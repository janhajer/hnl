#pragma once

#include <vector>
#include <string>
#include <map>
#include <functional>

#include <boost/filesystem/path.hpp>

namespace hnl {

std::string find_in_file_copy(std::vector<std::string>& lines, int pos, std::function<bool(std::vector<std::string>)> const& predicate);

std::vector<std::string> import_head(boost::filesystem::path const& path, int number);

std::vector<std::string> import_tail(boost::filesystem::path const& path, int number);

struct Meta {
    double mass = 0;
    double sigma = 0;
    std::map<int, std::map<int, double>> couplings;
};

}
