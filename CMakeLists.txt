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
add_library(hnl STATIC ResonanceWidths.cpp ThreeBodyWidth.cpp FormFactors.cpp MesonDecays.cpp BaryonDecays.cpp TauResonance.cpp branching_ratios.cpp)
target_include_directories(hnl PUBLIC "${PYTHIA8_INCLUDE_DIR}")
target_compile_definitions(hnl PUBLIC PYTHIA8_XMLDOC="${PYTHIA8_XMLDOC}")
target_link_libraries(hnl PUBLIC "${PYTHIA8_LIBRARY}" Boost::headers ${CMAKE_DL_LIBS})

# executables
foreach(name b-mesons resonance write-BRs write-lifetime)
    add_executable(${name} main-${name}.cpp)
    target_link_libraries(${name} hnl)
endforeach()

# hepmc3 (optional): -DHepMC3_DIR=<prefix>/share/HepMC3/cmake if not found
find_package(HepMC3 CONFIG QUIET)
if(HepMC3_FOUND)
    message(STATUS "HepMC3: ${HEPMC3_LIBRARIES}")
    add_executable(write-hepmc main-write-hepmc.cpp)
    target_include_directories(write-hepmc PRIVATE "${HEPMC3_INCLUDE_DIR}")
    target_link_libraries(write-hepmc hnl ${HEPMC3_LIBRARIES})
else()
    message(STATUS "HepMC3 not found: skipping write-hepmc")
endif()

# tests (Catch2 v3): installed, else downloaded at configure time; run with ctest, disable with -DBUILD_TESTING=OFF
include(CTest)
if(BUILD_TESTING)
    find_package(Catch2 3 QUIET)
    if(NOT Catch2_FOUND)
        include(FetchContent)
        FetchContent_Declare(Catch2 GIT_REPOSITORY https://github.com/catchorg/Catch2.git GIT_TAG v3.9.1 GIT_SHALLOW TRUE)
        FetchContent_MakeAvailable(Catch2)
        list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
    endif()
    add_executable(tests tests/form_factors.cpp tests/widths.cpp tests/tau.cpp)
    target_include_directories(tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    target_link_libraries(tests hnl Catch2::Catch2WithMain)
    include(Catch)
    catch_discover_tests(tests)
endif()
