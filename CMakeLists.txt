cmake_minimum_required (VERSION 2.8)

project(Neutrino)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -Ofast -Wall -Wextra -Wpedantic -Wfatal-errors")

find_package(ROOT)
find_package(Boost COMPONENTS system filesystem REQUIRED)

include_directories("${ROOT_INCLUDE_DIRS}")
include_directories("${Boost_INCLUDE_DIR}")

link_directories(~/local/lib/)
add_library(libDelphes SHARED IMPORTED)

add_executable (neutrino main.cpp)
target_link_libraries(neutrino Delphes ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
set_property(TARGET neutrino PROPERTY CXX_STANDARD 14)
