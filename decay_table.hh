#pragma once

#include <functional>
#include <vector>
#include <string>

namespace hnl {

std::vector<std::string> meson_decay_table(std::function<double (int id_heavy, int id_light)> const& coupling, double mass, int id, int meson_id);

std::vector<std::string> hnl_decay_table(std::function<double (int id_heavy, int id_light)> const& coupling, double mass, int id);

}
