cmake_minimum_required (VERSION 2.8)

project(Neutrino)

set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -Ofast -Wall -Wextra -Wpedantic -Wfatal-errors -Wno-language-extension-token -pthread")
# -g -ggdb

# message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(ROOT)

set(HepMC3_DIR ${PROJECT_SOURCE_DIR})
find_package(HepMC3)

find_package(Boost COMPONENTS system filesystem iostreams REQUIRED)

include_directories("${ROOT_INCLUDE_DIRS}")
include_directories("${Boost_INCLUDE_DIR}")

link_directories(~/local/lib/)
add_library(libDelphes SHARED IMPORTED)

set(libs ${Boost_LIBRARIES} ${ROOT_LIBRARIES})

set(delphes-libs ${libs} Delphes ${HepMC3_LIBRARIES})

set(hepmc-libs ${Boost_LIBRARIES} ${HepMC3_LIBRARIES})

set(hepmc-libs-plot ${libs} ${HepMC3_LIBRARIES})

# add_executable(CMS13 main-cms-13.cpp)
# target_link_libraries(CMS13 ${delphes-libs})
#
# add_executable(CMS14 main-cms-14.cpp)
# target_link_libraries(CMS14 ${delphes-libs})
#
# add_executable(ATLAS main-atlas.cpp)
# target_link_libraries(ATLAS ${delphes-libs})
#
# add_executable(LHCb main-lhcb.cpp)
# target_link_libraries(LHCb ${delphes-libs})
#
# add_executable(Distance main-distance.cpp)
# target_link_libraries(Distance ${delphes-libs})

add_executable(Test main-test.cpp)
target_link_libraries(Test ${libs})

# add_executable(Test2 main-ray.cpp)
# target_link_libraries(Test2 ${delphes-libs})
#
# add_executable(Wboson main-W.cpp)
# target_link_libraries(Wboson ${delphes-libs})
#
# add_executable(Leading main-lead.cpp)
# target_link_libraries(Leading ${delphes-libs})
#
# add_executable(HNL main-HNL.cpp)
# target_link_libraries(HNL ${delphes-libs})

# add_executable(Tau main-tau.cpp)
# target_link_libraries(Tau ${hepmc-libs-plot})

add_executable(CMS13-Pythia main-cms-13-pythia.cpp)
target_link_libraries(CMS13-Pythia ${hepmc-libs})

add_executable(CMS14-Pythia main-cms-14-pythia.cpp)
target_link_libraries(CMS14-Pythia ${hepmc-libs})

add_executable(ATLAS-Pythia main-atlas-pythia.cpp)
target_link_libraries(ATLAS-Pythia ${hepmc-libs})

add_executable(LHCb-Pythia main-lhcb-pythia.cpp)
target_link_libraries(LHCb-Pythia ${hepmc-libs})

# add_executable(SimpleTest maintest.cpp)
# target_link_libraries(SimpleTest ${libs})
#
# add_executable(VerticesATLAS main-vertices-atlas.cpp)
# target_link_libraries(VerticesATLAS ${hepmc-libs-plot})
#
# add_executable(VerticesCMS main-vertices-cms.cpp)
# target_link_libraries(VerticesCMS ${hepmc-libs-plot})
#
# add_executable(VerticesLHCb main-vertices-lhcb.cpp)
# target_link_libraries(VerticesLHCb ${hepmc-libs-plot})
#
# add_executable(Conservation main-conservation.cpp)
# target_link_libraries(Conservation ${hepmc-libs})

# macro(test name file)
# add_executable(${name} ${file})
# target_link_libraries(${name} ${hepmc-libs})
# add_test(NAME ${name} COMMAND ${name})
# endmacro()

# test(detector test.cpp)
