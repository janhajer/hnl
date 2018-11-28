cmake_minimum_required (VERSION 2.8)

project(Neutrino)

set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -Ofast -Wall -Wextra -Wpedantic -Wfatal-errors -Wno-language-extension-token -pthread")

# message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(ROOT)

set(HepMC3_DIR ${PROJECT_SOURCE_DIR})
find_package(HepMC3)

find_package(Boost COMPONENTS system filesystem iostreams REQUIRED)

include_directories("${ROOT_INCLUDE_DIRS}")
include_directories("${Boost_INCLUDE_DIR}")

link_directories(~/local/lib/)
add_library(libDelphes SHARED IMPORTED)

set(LIBRARIES ${Boost_LIBRARIES} ${ROOT_LIBRARIES} ${HepMC3_LIBRARIES} Delphes)

add_executable(CMS13 main-cms-13.cpp)
target_link_libraries(CMS13 ${LIBRARIES})

add_executable(CMS14 main-cms-14.cpp)
target_link_libraries(CMS14 ${LIBRARIES})

add_executable(ATLAS main-atlas.cpp)
target_link_libraries(ATLAS ${LIBRARIES})

add_executable(LHCb main-lhcb.cpp)
target_link_libraries(LHCb ${LIBRARIES})

add_executable(Distance main-distance.cpp)
target_link_libraries(Distance ${LIBRARIES})

add_executable(Test main-test.cpp)
target_link_libraries(Test ${LIBRARIES})

add_executable(Test2 main-ray.cpp)
target_link_libraries(Test2 ${LIBRARIES})

add_executable(Wboson main-W.cpp)
target_link_libraries(Wboson ${LIBRARIES})

add_executable(Leading main-lead.cpp)
target_link_libraries(Leading ${LIBRARIES})

add_executable(HNL main-HNL.cpp)
target_link_libraries(HNL ${LIBRARIES})

add_executable(CMS13-Pythia main-cms-13-pythia.cpp)
target_link_libraries(CMS13-Pythia ${LIBRARIES})

add_executable(CMS14-Pythia main-cms-14-pythia.cpp)
target_link_libraries(CMS14-Pythia ${LIBRARIES})

add_executable(ATLAS-Pythia main-atlas-pythia.cpp)
target_link_libraries(ATLAS-Pythia ${LIBRARIES})

add_executable(LHCb-Pythia main-lhcb-pythia.cpp)
target_link_libraries(LHCb-Pythia ${LIBRARIES})
