#pragma once

#include <iostream>

namespace hnl {

inline void print() {
    std::cout << std::endl;
}

template<typename Object, typename ... Arguments>
void print(Object const& object, Arguments ... arguments) {
    std::cout << std::boolalpha << std::scientific << object << ' ';
    print(arguments ...);
}

template<typename Container>
void print_line(Container const& container) {
    for (auto const& element : container) std::cout << element << ", ";
    std::cout << std::endl;
}

}
