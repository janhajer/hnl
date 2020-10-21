#pragma once

#include <iostream>

namespace hnl {

inline void print() {
    std::clog << std::endl;
}

template<typename Object, typename ... Arguments>
void print(Object const& object, Arguments ... arguments) {
    std::clog << std::boolalpha << std::scientific << object << ' ';
    print(arguments ...);
}

}
