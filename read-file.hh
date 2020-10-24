#pragma once

#include <vector>
#include <string>
#include <map>
#include <functional>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/range/iterator_range_core.hpp>

namespace hnl {

std::vector<std::string> split_line(std::string const& line);

std::string find_if(std::vector<std::string> const& lines, int pos, std::function<bool(std::vector<std::string>)> const& predicate);

std::vector<std::string> import_file(boost::filesystem::path const& path);

std::vector<std::string> import_head(boost::filesystem::path const& path, int number);

std::vector<std::string> import_tail(boost::filesystem::path const& path, int number);

using Couplings = std::map<int, std::map<int, double>>;

double totalvalue(Couplings const& couplings);

struct Meta {
    double mass = 0;
    double sigma = 0;
    Couplings couplings;
};

std::ostream& operator<<(std::ostream& stream, Meta const& meta);

double max(Couplings const& couplings);

using Result = std::map<double, std::map<double, double>>;

void save(Result const& result, std::string const& name = "result");

using Files = boost::iterator_range<boost::filesystem::directory_iterator>;

Files files(boost::filesystem::path const& path);

}
