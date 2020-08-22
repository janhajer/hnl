#pragma once

#include <map>
#include <tuple>

#include "math.hh"

namespace hnl {

using Result = std::map<int, std::map<std::tuple<int, int, int, int, int>, std::map<int, double>>>;

Result branching_ratio(Loop const& loop, double& mass_max, int source, int step);

void write_branching_ratios(int source);

void write_branching_ratios();

}

