cmake_minimum_required(VERSION 3.15)

project(HeavyNeutralLeptons CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_compile_options(-Wall -Wextra -Wpedantic -Wfatal-errors) # Release: -O3; no -Ofast, its -ffast-math breaks NaN and Inf handling
# add_compile_options(-march=native -mtune=native) # faster for local machine
# add_compile_options(-g -ggdb) # debugging

# pythia: -DPYTHIA8_DIR=<prefix>, else $PYTHIA8, else ~/Development/pythia
if(DEFINED ENV{PYTHIA8})
    set(pythia_default "$ENV{PYTHIA8}")
else()
    set(pythia_default "$ENV{HOME}/Development/pythia")
endif()
set(PYTHIA8_DIR "${pythia_default}" CACHE PATH "Pythia 8.3 installation prefix")
find_path(PYTHIA8_INCLUDE_DIR Pythia8/Pythia.h HINTS "${PYTHIA8_DIR}/include" NO_DEFAULT_PATH)
find_library(PYTHIA8_LIBRARY pythia8 HINTS "${PYTHIA8_DIR}/lib" NO_DEFAULT_PATH)
if(NOT PYTHIA8_INCLUDE_DIR OR NOT PYTHIA8_LIBRARY)
    message(FATAL_ERROR "Pythia 8.3 not found in ${PYTHIA8_DIR}; set -DPYTHIA8_DIR=<prefix>")
endif()
find_path(PYTHIA8_XMLDOC Index.xml HINTS "${PYTHIA8_DIR}/share/Pythia8/xmldoc" NO_DEFAULT_PATH)
if(NOT PYTHIA8_XMLDOC)
    message(FATAL_ERROR "Pythia data not found in ${PYTHIA8_DIR}/share/Pythia8/xmldoc")
endif()
message(STATUS "Pythia: ${PYTHIA8_LIBRARY}, data: ${PYTHIA8_XMLDOC}")

# boost: only header-only parts (range, optional, algorithm)
find_package(Boost REQUIRED CONFIG)

# lib
add_library(hnl STATIC ResonanceWidths.cpp ThreeBodyWidth.cpp MesonDecays.cpp BaryonDecays.cpp branching_ratios.cpp)
target_include_directories(hnl PUBLIC "${PYTHIA8_INCLUDE_DIR}")
target_compile_definitions(hnl PUBLIC PYTHIA8_XMLDOC="${PYTHIA8_XMLDOC}")
target_link_libraries(hnl PUBLIC "${PYTHIA8_LIBRARY}" Boost::headers ${CMAKE_DL_LIBS})

# executables
foreach(name b-mesons resonance write-BRs write-lifetime)
    add_executable(${name} main-${name}.cpp)
    target_link_libraries(${name} hnl)
endforeach()

# hepmc2 (optional): -DHEPMC2_DIR=<prefix>, needs Pythia8Plugins/HepMC2.h
set(HEPMC2_DIR "$ENV{HOME}/Development/HEPMC2" CACHE PATH "HepMC2 installation prefix")
find_path(HEPMC2_INCLUDE_DIR HepMC/GenEvent.h HINTS "${HEPMC2_DIR}" "${HEPMC2_DIR}/include")
find_library(HEPMC2_LIBRARY HepMC HINTS "${HEPMC2_DIR}/lib" "${HEPMC2_DIR}/build/lib")
if(HEPMC2_INCLUDE_DIR AND HEPMC2_LIBRARY)
    add_executable(write-hepmc main-write-hepmc.cpp)
    target_include_directories(write-hepmc PRIVATE "${HEPMC2_INCLUDE_DIR}")
    target_link_libraries(write-hepmc hnl "${HEPMC2_LIBRARY}")
else()
    message(STATUS "HepMC2 not found: skipping write-hepmc")
endif()
