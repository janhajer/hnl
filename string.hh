#pragma once

#include <string>
#include <sstream>

#include "io.hh"

namespace hnl {

using namespace std::string_literals;

inline std::string to_string(double value) {
    std::stringstream string_stream;
    string_stream << std::scientific << value;
    return string_stream.str();
}

inline double to_double(std::string const& string) {
    try {
        return std::stod(string);
    } catch (...) {
        print("The string:", string, ", is not a number");
        return 0.;
    }
}

inline auto filter(std::string string, std::string const& pattern) {
    auto position = string.find(pattern);
    while (position != std::string::npos) {
        string.erase(position, pattern.length());
        position = string.find(pattern);
    }
    return string;
}

template<typename Get, typename Condition, typename Second>
auto get_while_do(Get get_check, Condition condition, Second do_work) {
    auto check = get_check();
    while (condition(check)) {
        do_work(check);
        check = get_check();
    }
}

inline auto wfilter(std::string string, std::string const& pattern) {
    get_while_do([&]() {
        return string.find(pattern);
    }, [](auto position) {
        return position != std::string::npos;
    }, [&](auto position) {
        return string.erase(position, pattern.length());
    });
    return string;
}

}
