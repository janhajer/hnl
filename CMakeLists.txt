cmake_minimum_required(VERSION 3.4)

project(Neutrino)

include_directories("~/Development/pythia/include")
include_directories("~/Development/HEPMC2/")
LINK_DIRECTORIES("~/Development/pythia/lib")
LINK_DIRECTORIES("~/Development/HEPMC2/build/lib/")

set(graphic ON)

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -Ofast -Wall -Wextra -Wpedantic -Wfatal-errors -Wno-language-extension-token -pthread -fPIC -ldl")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -Wall -Wextra -Wpedantic -Wfatal-errors -Wno-language-extension-token -pthread -fPIC -ldl")



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

# hepmc
# set(HepMC3_DIR ${PROJECT_SOURCE_DIR})
# find_package(HepMC3)

# root
# find_package(ROOT)
# include_directories("${ROOT_INCLUDE_DIRS}")
# set(root-libs ${ROOT_LIBRARIES})

# delphes
# link_directories(~/local/lib/)
# add_library(libDelphes SHARED IMPORTED)


# set(libs ${boost-libs} ${ROOT_LIBRARIES})
#
# set(delphes-libs ${libs} Delphes ${HepMC3_LIBRARIES})
#

# boost
find_package(Boost COMPONENTS system filesystem iostreams REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
set(boost-libs ${Boost_LIBRARIES})

set(hepmc-libs ${boost-libs} ${HepMC3_LIBRARIES} ${cgal-libs})

set(plot-libs ${hepmc-libs} ${root-libs})

add_library(pythia ResonanceWidths.cpp ThreeBodyWidth.cpp Sigma.cpp)

add_executable(calculate-sigma main-calculate-sigma.cpp)
target_link_libraries(calculate-sigma pythia8 HepMC ${boost-libs} pythia)

add_executable(write-sigma-hepmc main-write-sigma-hepmc.cpp )
target_link_libraries(write-sigma-hepmc pythia8 HepMC ${boost-libs} pythia)

add_executable(write-sigma-hepmcs main-write-sigma-hepmcs.cpp )
target_link_libraries(write-sigma-hepmcs pythia8 HepMC ${boost-libs} pythia)

add_executable(write-BRs main-write-BRs.cpp )
target_link_libraries(write-BRs pythia8 HepMC ${boost-libs} pythia)

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

if(graphic)
  add_executable(Geometry geometry.cpp)
  target_link_libraries(Geometry HepMC ${hepmc-libs} ${qt-libs})
endif(graphic)
