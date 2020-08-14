#pragma once

#include <iostream>
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

namespace hnl
{

template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

auto arguments(int argc, char** argv) noexcept -> std::vector<std::string>
{
    return {argv + 1, argv + argc};
}

using namespace std::string_literals;

template<typename Integer>
auto irange(Integer integer) noexcept
{
    return boost::irange(0, integer);
}

template<typename Element>
auto make_vector(Element&& element) -> std::vector<Element>
{
    return {element};
}

template<typename Value>
auto to_string(std::vector<Value> const& vector) noexcept
{
    std::stringstream result;
    std::copy(vector.begin(), vector.end(), std::ostream_iterator<Value>(result, " "));
    return result.str();
}

// template <template<class...> class Container, typename Element, typename Function, typename Result = std::decay_t<std::result_of_t<Function&(Element const&)>>>
// std::vector<Result> transform(Container<Element> const& container, Function && function) noexcept {
//     std::vector<Result> result;
//     result.reserve(container.size());
//     boost::range::transform(container, std::back_inserter(result), function);
//     return result;
// }

template<typename Element>
auto& operator+=(std::vector<Element>& one, std::vector<Element> const& two) noexcept
{
    one.insert(one.end(), two.begin(), two.end());
    return one;
}

template<typename Element>
auto operator+(std::vector<Element> one, std::vector<Element> const& two) noexcept
{
//     auto copy = one;
    return one += two;
//     copy.insert(copy.end(), two.begin(), two.end());
//     return copy;
}

template<typename Element>
auto& operator+=(std::vector<Element>& one, Element const& two) noexcept
{
    one.emplace_back(two);
    return one;
}

template<typename Element>
auto size(std::vector<Element> const& vector) noexcept
{
    return vector.size();
}

template<typename Element>
auto operator!(std::vector<Element> const& container) noexcept
{
    return container.empty();
}

template<class Range, class Predicate>
auto has_at_least(Range const& range, int max, Predicate predicate) noexcept
{
    auto counter = 0;
    return boost::algorithm::any_of(range, [&](auto const & element) noexcept {
        return predicate(element) && ++counter >= max;
    });
}

template <template<class...> class Container, typename Element, typename Function, typename Result = std::decay_t<std::result_of_t<Function&(Element const&)>>>
          std::vector<Result> transform(Container<Element> const& container, Function && function) noexcept
{
    std::vector<Result> result;
    result.reserve(size(container));
    boost::range::transform(container, std::back_inserter(result), function);
    return result;
}

template<template<class...> class Container, typename Element, typename Function>
auto copy_if(Container<Element> const& container, Function function) noexcept
{
    std::vector<Element> result(size(container));
    auto iterator = boost::algorithm::copy_if(container, std::begin(result), function);
    result.erase(iterator, result.end());
    return result;
}

template<typename Element, typename Predicate>
auto find_erase(std::vector<Element>& container, Predicate predicate) noexcept -> boost::optional<Element>
{
    auto found = boost::range::find_if(container, predicate);
    if (found == container.end()) return boost::none;
    auto element = *found;
    container.erase(found);
    return element;
}

template<typename Object>
auto sqr(Object const& object) noexcept
{
    return object * object;
}

template<typename Object>
auto cube(Object const& object) noexcept
{
    return object * object * object;
}

template<typename Object>
auto norm(Object const& one) noexcept
{
    return sqrt(sqr(one));
}

template<typename Key_, typename Value_>
auto& operator<<(std::ostream& stream, std::pair<Key_, Value_> const& pair) noexcept
{
    return stream << '(' << pair.first << ", " << pair.second << ')';
}

template<typename Element, template <typename, typename = std::allocator<Element>> class Container>
auto & operator<<(std::ostream& stream, Container<Element> const& container) noexcept
{
    for (auto const& element : boost::adaptors::index(container)) stream << '\n' << element.index() << ": " << element.value();
    return stream;
}

void print() noexcept
{
    std::cout << std::endl;
}

template<typename Object, typename ... Arguments>
void print(Object const& object, Arguments ... arguments) noexcept
{
    std::cout << std::boolalpha << std::scientific << object << ' ';
    print(arguments ...);
}

template<typename Container>
void print_line(Container const& container) noexcept
{
    for (auto const& element : container) std::cout << element << ", ";
    std::cout << std::endl;
}

auto filter(std::string string, std::string const& pattern) noexcept
{
    auto position = string.find(pattern);
    while (position != std::string::npos) {
        string.erase(position, pattern.length());
        position = string.find(pattern);
    }
    return string;
}

template <typename Enumeration>
constexpr typename std::underlying_type_t<Enumeration> to_underlying(Enumeration enumeration) noexcept
{
    return static_cast<typename std::underlying_type_t<Enumeration>>(enumeration);
}

template <typename Enumeration>
std::string to_string(Enumeration enumeration) noexcept
{
    return std::to_string(to_underlying(enumeration));
}

template <typename Element>
auto distinct_pairs(std::vector<Element> const& vector) noexcept
{
    std::vector<std::pair<Element, Element>> pairs;
    for (auto i = 0u; i < vector.size(); ++i) for (auto j = i + 1u; j < vector.size(); ++j) pairs.emplace_back(vector.at(i), vector.at(j));
    return pairs;
}

template<typename Get, typename Condition, typename Second>
auto get_while_do(Get get_check, Condition condition, Second do_work)
{
    auto check = get_check();
    while (condition(check)) {
        do_work(check);
        check = get_check();
    }
}

auto wfilter(std::string string, std::string const& pattern) noexcept
{
    get_while_do([&]() noexcept {
        return string.find(pattern);
    }, [](auto position) noexcept {
        return position != std::string::npos;
    }, [&](auto position) noexcept {
        return string.erase(position, pattern.length());
    });
    return string;
}

template<typename Container, typename Predicate>
auto& keep_if(Container& container, Predicate&&  predicate) noexcept
{
    return boost::range::remove_erase_if(container, [predicate](auto const & element) noexcept {
        return !predicate(element);
    });
}

template<typename Base>
Base pow(Base base, int exp)
{
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
