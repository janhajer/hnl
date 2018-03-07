cmake_minimum_required (VERSION 2.8)

project(Neutrino)

find_package(ROOT)
find_package(Boost)
include_directories("${ROOT_INCLUDE_DIRS}")

add_library(libDelphes SHARED IMPORTED)

add_executable (neutrino main.cpp)
target_link_libraries(neutrino Delphes ${ROOT_LIBRARIES})
set_property(TARGET neutrino PROPERTY CXX_STANDARD 14)
