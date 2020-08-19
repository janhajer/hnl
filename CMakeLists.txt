cmake_minimum_required(VERSION 3.4)

project(HeavyNeutralLeptons)

include_directories("~/Development/pythia/include")
include_directories("~/Development/HEPMC2/")
LINK_DIRECTORIES("~/Development/pythia/lib")
LINK_DIRECTORIES("~/Development/HEPMC2/build/lib/")

set(graphic OFF)

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -Ofast -Wall -Wextra -Wpedantic -Wfatal-errors -Wno-language-extension-token -pthread -fPIC -ldl")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -Wall -Wextra -Wpedantic -Wfatal-errors -Wno-language-extension-token -pthread -fPIC -ldl -lz")

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -ggdb")

# boost
find_package(Boost COMPONENTS system filesystem iostreams REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
set(boost-libs ${Boost_LIBRARIES})

# cgal
if(graphic)
  find_package(CGAL REQUIRED COMPONENTS Core Qt5)
  set(cgal-libs ${CGAL_LIBRARIES})
elseif(graphic)
  set(cgal-libs )
endif(graphic)

# qt
if(graphic)
  find_package(Qt5 COMPONENTS Widgets Gui Core OpenGL Xml REQUIRED)
  include_directories(${Qt5Widgets_INCLUDE_DIRS})
  include_directories(${Qt5Xml_INCLUDE_DIRS})
  include_directories(${Qt5OpenGL_INCLUDE_DIRS})
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTOUIC_SEARCH_PATHS /usr/include/CGAL/Qt/resources)
  set(qt-libs ${Qt5Widgets_LIBRARIES})
elseif(graphic)
  set(qt-libs )
endif(graphic)

# boost
find_package(Boost COMPONENTS system filesystem iostreams REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
set(boost-libs ${Boost_LIBRARIES})

# lib
add_library(pythia ResonanceWidths.cpp ThreeBodyWidth.cpp Sigma.cpp)

# executables
add_executable(resonance main-resonance.cpp)
target_link_libraries(resonance pythia8 HepMC ${boost-libs} pythia)

add_executable(calculate-sigma main-calculate-sigma.cpp)
target_link_libraries(calculate-sigma pythia8 HepMC ${boost-libs} pythia)

add_executable(write-sigma-hepmc main-write-sigma-hepmc.cpp )
target_link_libraries(write-sigma-hepmc pythia8 HepMC ${boost-libs} pythia)

add_executable(write-sigma-hepmcs main-write-sigma-hepmcs.cpp )
target_link_libraries(write-sigma-hepmcs pythia8 HepMC ${boost-libs} pythia)

add_executable(write-BRs main-write-BRs.cpp )
target_link_libraries(write-BRs pythia8 ${boost-libs} pythia)

add_executable(write-hepmc main-write-hepmc.cpp )
target_link_libraries(write-hepmc pythia8 HepMC ${boost-libs} pythia)

add_executable(write-hepmcs main-write-hepmcs.cpp )
target_link_libraries(write-hepmcs pythia8 HepMC ${boost-libs} pythia)

add_executable(read-hepmc main-read-hepmc.cpp )
target_link_libraries(read-hepmc pythia8 HepMC ${boost-libs} pythia)

add_executable(scan-hepmc main-scan-hepmc.cpp )
target_link_libraries(scan-hepmc pythia8 HepMC ${boost-libs} pythia)

add_executable(scan-hepmcs main-scan-hepmcs.cpp )
target_link_libraries(scan-hepmcs pythia8 HepMC ${boost-libs} pythia)

add_executable(scan-lhe main-scan-lhe.cpp )
target_link_libraries(scan-lhe pythia8 HepMC ${boost-libs} pythia)

add_executable(scan-lhes main-scan-lhes.cpp )
target_link_libraries(scan-lhes pythia8 HepMC ${boost-libs} pythia)

add_executable(read-lhe main-read-lhe.cpp )
target_link_libraries(read-lhe pythia8 ${boost-libs} pythia)

if(graphic)
  add_executable(Geometry geometry.cpp)
  target_link_libraries(Geometry HepMC ${hepmc-libs} ${qt-libs})
endif(graphic)
