#pragma once

#include <functional>
#include <vector>
#include <string>
#include <map>

namespace hnl {

using BranchingRatios = std::map<int, std::map<std::array<int, 5>, std::map<int, double>>>;

std::vector<std::string> meson_decay_table(std::function<double (int id_heavy, int id_light)> const& coupling, double mass, int id, int meson_id);

std::vector<std::string> hnl_decay_table(std::function<double (int id_heavy, int id_light)> const& coupling, double mass, int id);

BranchingRatios hnl_branching_ratios(std::function<double (int id_heavy, int id_light)> const& coupling, double mass, int id, int step);

}
