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

// class ThreadPool
// {
//     public:
//     template<typename Index, typename Function>
// static void ParallelFor(Index start, Index end, Function function) noexcept {
//     // Estimate number of threads in the pool
//     static auto const nb_threads_hint = std::thread::hardware_concurrency();
//     static auto const nb_threads = nb_threads_hint == 0u ? 8u : nb_threads_hint;
//     // Size of a slice for the range functions
//     Index n = end - start + 1;
//     Index slice = static_cast<Index>(std::round(n / static_cast<double>(nb_threads)));
//     slice = std::max(slice, Index(1));
//     // [Helper] Inner loop
//     auto launchRange = [&function](int begin, int end) noexcept {
//         for (auto counter = begin; counter < end; ++counter) function(counter);
//     };
//     // Create pool and launch jobs
//     std::vector<std::thread> pool;
//     pool.reserve(nb_threads);
//     Index i1 = start;
//     Index i2 = std::min(start + slice, end);
//     for (auto i = 0u; i + 1 < nb_threads && i1 < end; ++i)
//     {
//         pool.emplace_back(launchRange, i1, i2);
//         i1 = i2;
//         i2 = std::min(i2 + slice, end);
//     }
//     if (i1 < end) pool.emplace_back(launchRange, i1, end);
//     for (std::thread& t : pool) if (t.joinable()) t.join(); // Wait for jobs to finish
// }
//
// // Serial version for easy comparison
// template<typename Index, typename Function>
// static void SequentialFor(Index start, Index end, Function function) noexcept {
//     for (auto i = start; i < end; i++) function(i);
// }
//     };




// auto example()
// {
//     std::mutex mutex;
//     ThreadPool::ParallelFor(0, 16, [&](int i) {
//         std::lock_guard<std::mutex> lock(mutex);
//         std::cout << i << std::endl;
//     });
// }

// template<typename Range, typename Function>
// auto parallel_for(Range const& range, Function function)
// {
// //     std::mutex mutex;
//     std::for_each(std::execution::par, range.begin(),range.end(),[](auto const& element) {
//         auto res += function(element)
// //         std::lock_guard<std::mutex> lock(mutex);
//
//     });
// }


// template <typename Iterator, typename Function>
// auto parallel_sum(Iterator begin, Iterator end, Function function) {
//     auto length = end - begin;
//     if(length < 1000) return std::accumulate(begin, end, 0);
//     Iterator center = begin + length/2;
//     auto async = std::async(std::launch::async, parallel_sum<Iterator, Function>, center, end, function);
//     auto sum = parallel_sum(begin, center, function);
//     return sum + async.get();
// }
//
// template <typename Range, typename Function>
// auto parallel_sum(Range const& range, Function function) {
//     return parallel_sum(range.begin(), range.end(), function);
// }

}
