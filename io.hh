#pragma once

#include <iostream>

namespace hnl {

inline void print() noexcept {
    std::cout << std::endl;
}

template<typename Object, typename ... Arguments>
void print(Object const& object, Arguments ... arguments) noexcept {
    std::cout << std::boolalpha << std::scientific << object << ' ';
                                 print(arguments ...);
}

template<typename Container>
void print_line(Container const& container) noexcept {
    for (auto const& element : container) std::cout << element << ", ";
            std::cout << std::endl;
}

}
