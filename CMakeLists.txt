cmake_minimum_required(VERSION 3.4)

project(HeavyNeutralLeptons)

# hepmc
include_directories("~/Development/HEPMC2/")
LINK_DIRECTORIES("~/Development/HEPMC2/build/lib/")

# pythia
include_directories("~/Development/pythia/include")
LINK_DIRECTORIES("~/Development/pythia/lib")

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_STANDARD 17) # we want this
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17") # but sometimes only this works
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Ofast -Wall -Wextra -Wpedantic -Wfatal-errors -Wno-language-extension-token")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native") # faster for local machine
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -ggdb") # debugging

# boost
find_package(Boost COMPONENTS system filesystem iostreams REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
set(boost-libs ${Boost_LIBRARIES})

# lib
add_library(hnl ResonanceWidths.cpp ThreeBodyWidth.cpp decay_table.cpp branching_ratios.cpp)

# libs
set(hnl-libs hnl pythia8 ${boost-libs} ${CMAKE_DL_LIBS} z)
set(hnl-libs-hep ${hnl-libs} HepMC)

# executables
add_executable(resonance main-resonance.cpp)
target_link_libraries(resonance ${hnl-libs})

add_executable(write-BRs main-write-BRs.cpp)
target_link_libraries(write-BRs ${hnl-libs})

add_executable(write-lifetime main-write-lifetime.cpp)
target_link_libraries(write-lifetime ${hnl-libs-hep})

add_executable(write-hepmc main-write-hepmc.cpp)
target_link_libraries(write-hepmc ${hnl-libs-hep})
