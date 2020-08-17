#pragma once

#include <vector>

#include "container.hh"
#include "io.hh"

namespace hnl {

const int heavy_neutrino = 9900012;

inline std::vector<int> heavy_neutral_leptons() {
    return {9900012, 9900014, 9900016};
}

inline bool is_heavy_neutral_lepton(int id) {
    return id == 9900012 || id == 9900014 || id == 9900016;
}

inline std::vector<int> light_neutrinos() {
    return {12, 14, 16};
}

inline std::vector<int> neutral_leptons() {
    return {12, 14, 16};
}

inline std::vector<int> charge_leptons() {
    return {11, 13, 15};
}

inline std::vector<int> leptons() {
    return neutral_leptons() + charge_leptons();
}

inline std::vector<int> down_type() {
    return {1, 3, 5};
}

inline std::vector<int> up_type() {
    return {2, 4, 6};
}

inline std::vector<int> quarks() {
    return down_type() + up_type();
}

inline std::vector<int> fermions() {
    return quarks() + leptons();
}

inline bool is_lepton(int id) {
    return id > 10 && id < 20;
}

inline bool is_charge_lepton(int id) {
    return is_lepton(id) && id % 2 == 1;
}

inline bool is_light_neutrino(int id) {
    return is_lepton(id) && id % 2 == 0;
}

inline bool is_quark(int id) {
    return id < 10;
}

inline bool is_up_type(int id) {
    return is_quark(id) && id % 2 == 0;
}

inline bool is_down_type(int id) {
    return is_quark(id) && id % 2 == 1;
}

inline bool is_meson(int id) {
    return id > 100 && id < 1000;
}

inline bool is_vector(int id) {
    return id == 113 || id == 213 || id == 313 || id == 323 || id == 413 || id == 423 || id == 433;
}

inline bool is_B_c_star(int id) {
    return id == 543;
}

inline bool is_B_c(int id) {
    return id == 541;
}

inline bool is_B_s_star(int id) {
    return id == 533;
}

inline bool is_B_s(int id) {
    return id == 531;
}

inline bool is_B_star(int id) {
    return id == 513 || id == 523;
}

inline bool is_B(int id) {
    return id == 511 || id == 521;
}

inline bool is_D_s_star(int id) {
    return id == 433;
}

inline bool is_D_s(int id) {
    return id == 431;
}

inline bool is_D_star(int id) {
    return id == 413 || id == 423;
}

inline bool is_D(int id) {
    return id == 411 || id == 421;
}

inline bool is_K_star(int id) {
    return id == 313 || id == 323 || id == 311 || id == 321;
}

inline bool is_K(int id) {
    return id == 130 || id == 310 || id == 311 || id == 321;
}

inline bool is_rho(int id) {
    return id == 113 || id == 213;
}

inline bool is_pi(int id) {
    return id == 111 || id == 211;
}

inline bool is_eta(int id) {
    return id == 221 || id == 331;
}

inline bool is_neutral_Kaon(int id) {
    return id == 130 || id == 310 || id == 311;
}

auto neutrino_coupling = [](double factor) {
    return [factor](int id_heavy, int id_light) -> double {
        if (!is_heavy_neutral_lepton(id_heavy)) {
            print(id_heavy, "is not a heavy neutrino");
            return 0;
        }
        if (!is_light_neutrino(id_light)) {
            print(id_light, "is not a light neutrino");
            return 0;
        }
        auto id = id_heavy - 9900000 - id_light;
        return id == 0 ? factor : 0;
        return id == 0 || std::abs(id) == 2 || std::abs(id) == 4 ? factor : 0;
    };
};

}
