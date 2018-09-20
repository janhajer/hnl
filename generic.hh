#pragma once

#include <iostream>
#include <vector>

#include <boost/units/cmath.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <boost/algorithm/cxx11/copy_if.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/range/size.hpp>

using namespace std::string_literals;

template<typename Integer>
auto irange(Integer integer) noexcept {
    return boost::irange(0, integer);
}

template<typename Element>
auto operator+(std::vector<Element> const& one, std::vector<Element> const& two) noexcept
{
    auto copy = one;
    copy.insert(copy.end(), two.begin(), two.end());
    return copy;
}

template<typename Element>
auto& operator+=(std::vector<Element>& one, std::vector<Element> const& two) noexcept {
    one.insert(one.end(), two.begin(), two.end());
    return one;
}

template<typename Element>
auto size(std::vector<Element> const& vector){
    return vector.size();
}

template<typename> class TTreeReaderArray;

template<typename Element>
auto size(TTreeReaderArray<Element> const& vector){
    return vector.GetSize();
}

template <template<class...> class Container, typename Element, typename Function, typename Result = std::decay_t<std::result_of_t<Function&(Element const&)>>>
std::vector<Result> transform(Container<Element> const& container, Function && function) noexcept {
    std::vector<Result> result;
    result.reserve(size(container));
    boost::range::transform(container, std::back_inserter(result), function);
    return result;
}

template<template<class...> class Container, typename Element, typename Function>
auto copy_if(Container<Element> const& container, Function function) noexcept {
    std::vector<Element> result(size(container));
    auto iterator = boost::algorithm::copy_if(container, std::begin(result), function);
    result.erase(iterator, result.end());
    return result;
}

template<typename Element, typename Predicate>
auto find_erase(std::vector<Element>& container, Predicate predicate) noexcept -> boost::optional<Element> {
    auto found = boost::range::find_if(container, predicate);
    if (found == container.end()) return boost::none;
    auto element = *found;
    container.erase(found);
    return element;
}

template<typename Object>
auto sqr(Object const& object) noexcept {
    return object * object;
}

template<typename Object>
auto length(Object const& one) noexcept {
    return sqrt(sqr(one));
}

template<typename Key_, typename Value_>
auto& operator<<(std::ostream& stream, std::pair<Key_, Value_> const& pair) noexcept {
    return stream << '(' << pair.first << ", " << pair.second << ')';
}

template<typename Element, template <typename, typename = std::allocator<Element>> class Container>
auto & operator<<(std::ostream& stream, Container<Element> const& container) noexcept {
    for (auto const& element : boost::adaptors::index(container)) stream << '\n' << element.index() << ": " << element.value();
    return stream;
}

void print() noexcept {
    std::cout << std::endl;
}

template<typename Object, typename ... Arguments>
void print(Object const& object, Arguments ... arguments) noexcept {
    std::cout << std::boolalpha << object << ' ';
    print(arguments ...);
}

template<typename Container>
void print_line(Container const& container) noexcept {
    for (auto const& element : container) std::cout << element << ", ";
    std::cout << std::endl;
}
