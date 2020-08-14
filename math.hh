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

template<typename Base>
Base pow(Base base, int exp) {
    assert(exp >= 0);
    Base result = 1;
    for (;;) {
        if (exp & 1) result *= base;
        exp >>= 1;
        if (!exp) break;
        base *= base;
    }
    return result;
}

inline auto lin_scale(double min, double max, int step, int steps) {
    return min + (max - min) * step / steps;
}

inline auto log_scale(double min, double max, int step, int steps) {
    return std::pow(10, lin_scale(std::log10(min), std::log10(max), step, steps));
}

}
