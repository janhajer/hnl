cmake_minimum_required (VERSION 2.8)

project(Neutrino)

find_package(ROOT)
find_package(Boost)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_library(libDelphes SHARED IMPORTED)

add_executable (neutrino main.cpp)
target_link_libraries(neutrino Delphes ${ROOT_LIBRARIES})
