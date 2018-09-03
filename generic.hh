#pragma once

#include <iostream>
#include <fstream>
#include <random>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/size.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <boost/optional/optional_io.hpp>

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

#include "classes/DelphesClasses.h"

#include "alphanum.hpp"
#include "detector.hh"

using namespace std::string_literals;

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

template<typename Particle>
auto transverse_distance(Particle const& particle) noexcept {
    return std::sqrt(sqr(particle.X) + sqr(particle.Y));
}

std::ostream& operator<<(std::ostream& stream, boost::filesystem::path const& path) noexcept {
    return stream << path.string();
}

std::ostream& operator<<(std::ostream& stream, GenParticle const& particle) noexcept {
    return stream << "ID: " << particle.PID << ", D0: " << transverse_distance(particle) << " mm" << " Status: " << particle.Status;
}

std::ostream& operator<<(std::ostream& stream, Muon const& muon) noexcept {
    return stream << "(" << muon.PT << ", " << muon.Eta << ")";
}

std::ostream& operator<<(std::ostream& stream, Electron const& muon) noexcept {
    return stream << "(" << muon.PT << ", " << muon.Eta << ")";
}

std::ostream& operator<<(std::ostream& stream, Jet const& muon) noexcept {
    return stream << "(" << muon.PT << ", " << muon.Eta << ")";
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
