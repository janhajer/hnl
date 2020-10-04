#pragma once

#include <vector>

#include <boost/range/irange.hpp>
#include <boost/range/algorithm/transform.hpp>

#include "math.hh"

namespace hnl {

template<typename Integer>
auto irange(Integer integer) {
    return boost::irange(static_cast<Integer>(0), integer);
}


template <template<class...> class Container, typename Element, typename Function, typename Result = std::decay_t<std::result_of_t<Function&(Element const&)>>>
std::vector<Result> transform(Container<Element> const& container, Function && function) {
    std::vector<Result> result;
    result.reserve(size(container));
    boost::range::transform(container, std::back_inserter(result), function);
    return result;
}

inline auto log_range(double min, double max, int steps) {
    return transform(irange(steps + 1), [&](auto step) {
        return log_step(min, max, step, steps);
    });
}

}

