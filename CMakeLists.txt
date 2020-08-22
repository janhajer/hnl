cmake_minimum_required(VERSION 3.4)

project(HeavyNeutralLeptons)

# use QT and cgal libs?
set(graphic ON)

# HepMC
include_directories("~/Development/HEPMC2/")
LINK_DIRECTORIES("~/Development/HEPMC2/build/lib/")

# pythia
include_directories("~/Development/pythia/include")
LINK_DIRECTORIES("~/Development/pythia/lib")

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_STANDARD 17) # we want this
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17") # but sometimes only this works
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -Wall -Wextra -Wpedantic -Wfatal-errors -Wno-language-extension-token -pthread -fPIC")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native") # faster for local machine
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -ggdb") # debugging

# cgal
if(graphic)
  find_package(CGAL REQUIRED COMPONENTS Core Qt5)
  set(cgal-libs ${CGAL_LIBRARIES})
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
endif(graphic)

# boost
find_package(Boost COMPONENTS system filesystem iostreams REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
set(boost-libs ${Boost_LIBRARIES})

# lib
add_library(hnl ResonanceWidths.cpp ThreeBodyWidth.cpp Sigma.cpp decay_table.cpp read-file.cpp branching_ratios.cpp)

set(hnl-libs hnl pythia8 ${boost-libs} ${CMAKE_DL_LIBS} z)
set(hnl-libs-hep ${hnl-libs} HepMC)

# executables
add_executable(resonance main-resonance.cpp)
target_link_libraries(resonance ${hnl-libs})

add_executable(calculate-sigma main-calculate-sigma.cpp)
target_link_libraries(calculate-sigma ${hnl-libs-hep})

add_executable(write-sigma-hepmc main-write-sigma-hepmc.cpp)
target_link_libraries(write-sigma-hepmc ${hnl-libs-hep})

add_executable(write-sigma-hepmcs main-write-sigma-hepmcs.cpp)
target_link_libraries(write-sigma-hepmcs ${hnl-libs-hep})

add_executable(write-BRs main-write-BRs.cpp)
target_link_libraries(write-BRs ${hnl-libs})

add_executable(write-hepmc main-write-hepmc.cpp)
target_link_libraries(write-hepmc ${hnl-libs-hep})

add_executable(write-hepmcs main-write-hepmcs.cpp)
target_link_libraries(write-hepmcs ${hnl-libs-hep})

add_executable(read-hepmc main-read-hepmc.cpp)
target_link_libraries(read-hepmc ${hnl-libs-hep})

add_executable(scan-hepmc main-scan-hepmc.cpp)
target_link_libraries(scan-hepmc ${hnl-libs-hep})

add_executable(scan-hepmcs main-scan-hepmcs.cpp)
target_link_libraries(scan-hepmcs ${hnl-libs-hep})

add_executable(scan-lhe main-scan-lhe.cpp)
target_link_libraries(scan-lhe ${hnl-libs})

add_executable(scan-lhes main-scan-lhes.cpp)
target_link_libraries(scan-lhes ${hnl-libs})

add_executable(read-lhe main-read-lhe.cpp)
target_link_libraries(read-lhe ${hnl-libs})

if(graphic)
  add_executable(Geometry main-geometry.cpp)
  target_link_libraries(Geometry ${qt-libs})

  add_executable(draw main-draw.cpp)
  target_link_libraries(draw ${qt-libs} ${cgal-libs})
endif(graphic)
