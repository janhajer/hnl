#pragma once

#include <cmath>
#include <cassert>

namespace hnl {

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename Object>
auto sqr(Object const& object) {
    return object * object;
}

template<typename Object>
auto cube(Object const& object) {
    return object * object * object;
}

template<typename Object>
auto norm(Object const& one) {
    return sqrt(sqr(one));
}

inline auto linear_value(double min, double max, int step, int steps) {
    return min + (max - min) * step / steps;
}

inline auto log_value(double min, double max, int step, int steps) {
    return std::pow(10, linear_value(std::log10(min), std::log10(max), step, steps));
}

inline auto linear_step(double min, double max, double point, int steps){
    return steps * (point - min) / (max - min);
}

inline auto log_step(double min, double max, double value, int steps) {
    return linear_step(std::log10(min), std::log10(max), std::log10(value), steps);
}


struct Loop {
    Loop(double min, int steps_) : m_min(min), steps(steps_) {}
    double m_min;
    int steps;
    double mass(double max, int step) const {
        return log_value(m_min, max, step, steps);
    }
};

}
