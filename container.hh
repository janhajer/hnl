#pragma once

#include <map>
#include <vector>
#include <sstream>

#include <boost/optional/optional_io.hpp>
#include <boost/algorithm/cxx11/copy_if.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/range/size.hpp>

#include "math.hh"

namespace hnl {

template<typename Integer>
auto irange(Integer integer) {
    return boost::irange(0, integer);
}

template<typename Element>
auto make_vector(Element&& element) -> std::vector<Element> {
    return {element};
}

template<typename Value>
auto to_string(std::vector<Value> const& vector) {
    std::stringstream result;
    std::copy(vector.begin(), vector.end(), std::ostream_iterator<Value>(result, " "));
    return result.str();
}

// template <template<class...> class Container, typename Element, typename Function, typename Result = std::decay_t<std::result_of_t<Function&(Element const&)>>>
// std::vector<Result> transform(Container<Element> const& container, Function && function) {
//     std::vector<Result> result;
//     result.reserve(container.size());
//     boost::range::transform(container, std::back_inserter(result), function);
//     return result;
// }

template<typename Element>
auto& operator+=(std::vector<Element>& one, std::vector<Element> const& two) {
    one.insert(one.end(), two.begin(), two.end());
    return one;
}

template<typename Element>
auto operator+(std::vector<Element> one, std::vector<Element> const& two) {
//     auto copy = one;
    return one += two;
//     copy.insert(copy.end(), two.begin(), two.end());
//     return copy;
}

template<typename Element>
auto& operator+=(std::vector<Element>& one, Element const& two) {
    one.emplace_back(two);
    return one;
}

template<typename Element>
auto size(std::vector<Element> const& vector) {
    return vector.size();
}

template<typename Element>
auto operator!(std::vector<Element> const& container) {
    return container.empty();
}

template<class Range, class Predicate>
auto has_at_least(Range const& range, int max, Predicate predicate) {
    auto counter = 0;
    return boost::algorithm::any_of(range, [&](auto const & element) {
        return predicate(element) && ++counter >= max;
    });
}

template <template<class...> class Container, typename Element, typename Function, typename Result = std::decay_t<std::result_of_t<Function&(Element const&)>>>
std::vector<Result> transform(Container<Element> const& container, Function && function) {
    std::vector<Result> result;
    result.reserve(size(container));
    boost::range::transform(container, std::back_inserter(result), function);
    return result;
}

template<template<class...> class Container, typename Element, typename Function>
auto copy_if(Container<Element> const& container, Function function) {
    std::vector<Element> result(size(container));
    auto iterator = boost::algorithm::copy_if(container, std::begin(result), function);
    result.erase(iterator, result.end());
    return result;
}

template<typename Element, typename Predicate>
auto find_erase(std::vector<Element>& container, Predicate predicate) -> boost::optional<Element> {
    auto found = boost::range::find_if(container, predicate);
    if (found == container.end()) return boost::none;
    auto element = *found;
    container.erase(found);
    return element;
}

template<typename Key_, typename Value_>
auto& operator<<(std::ostream& stream, std::pair<Key_, Value_> const& pair) {
    return stream << '(' << pair.first << ", " << pair.second << ')';
}

template<typename Key, typename Value>
auto & operator<<(std::ostream& stream, std::map<Key, Value> const& map) {
    for (auto const& pair : map) stream << '\n' << pair.first << ": " << pair.second;
    return stream;
}

template<typename Element, template <typename, typename = std::allocator<Element>> class Container>
auto & operator<<(std::ostream& stream, Container<Element> const& container) {
    for (auto const& element : boost::adaptors::index(container)) stream << '\n' << element.index() << ": " << element.value();
    return stream;
}

template<typename Element, size_t size>
auto& operator<<(std::ostream& stream, std::array<Element, size> const& container) {
    for (auto const& element : boost::adaptors::index(container)) stream << '\n' << element.index() << ": " << element.value();
    return stream;
}

template<typename Container, typename Predicate>
auto& keep_if(Container& container, Predicate&&  predicate) {
    return boost::range::remove_erase_if(container, [predicate](auto const & element) {
        return !predicate(element);
    });
}

template<typename First, typename Second>
std::map<First, Second>& operator+=(std::map<First, Second>& left, std::map<First, Second> const& right) {
    for (auto const& pair : right) left[pair.first] += pair.second;
    return left;
}

inline auto log_range(double min, double max, int steps) {
    return transform(irange(steps + 1), [&](auto step) {
        return log_scale(min, max, step, steps);
    });
}

}

