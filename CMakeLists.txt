cmake_minimum_required(VERSION 3.4)

project(Neutrino)

include_directories("~/Development/pythia/include")
include_directories("~/Development/HEPMC2/")
LINK_DIRECTORIES("~/Development/pythia/lib")
LINK_DIRECTORIES("~/Development/HEPMC2/build/lib/")

set(graphic OFF)

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

# message(${hepmc-libs})
#
set(plot-libs ${hepmc-libs} ${root-libs})

# add_library(pythia-test ResonanceWidths.cpp ThreeBodyWidth.cpp)

add_executable(MAPP main-pythia.cpp ResonanceWidths.cpp ThreeBodyWidth.cpp)
target_link_libraries(MAPP pythia8 HepMC ${boost-libs})

# add_executable(Mapp main-mapp-pythia.cpp)
# target_link_libraries(Mapp ${hepmc-libs})
#
# add_executable(MappPlot main-mapp-plot.cpp)
# target_link_libraries(MappPlot ${plot-libs})

if(graphic)
  add_executable(Geometry geometry.cpp)
  target_link_libraries(Geometry HepMC ${hepmc-libs} ${qt-libs})
endif(graphic)

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

# add_executable(Test main-test.cpp)
# target_link_libraries(Test ${libs})

# add_executable(Test-Eff main-eff.cpp)
# target_link_libraries(Test-Eff ${plot-libs})

# add_executable(Test2 main-ray.cpp)
# target_link_libraries(Test2 ${delphes-libs})
#
# add_executable(Wboson main-W.cpp)
# target_link_libraries(Wboson ${plot-libs})
#
# add_executable(PhaseSpace main-PS.cpp)
# target_link_libraries(PhaseSpace ${plot-libs})
#
# add_executable(Leading main-lead.cpp)
# target_link_libraries(Leading ${delphes-libs})
#
# add_executable(HNL main-HNL.cpp)
# target_link_libraries(HNL ${delphes-libs})

# add_executable(Tau main-tau.cpp)
# target_link_libraries(Tau ${plot-libs})

# add_executable(CMS13-Pythia main-cms-13-pythia.cpp)
# target_link_libraries(CMS13-Pythia ${hepmc-libs})
#
# add_executable(CMS14-Pythia main-cms-14-pythia.cpp)
# target_link_libraries(CMS14-Pythia ${hepmc-libs})
#
# add_executable(ATLAS-Pythia main-atlas-pythia.cpp)
# target_link_libraries(ATLAS-Pythia ${hepmc-libs})
#
# add_executable(LHCb-Pythia main-lhcb-pythia.cpp)
# target_link_libraries(LHCb-Pythia ${hepmc-libs})
#
# add_executable(MAPP-Pythia main-mapp-pythia.cpp)
# target_link_libraries(MAPP-Pythia ${hepmc-libs})

# add_executable(SimpleTest maintest.cpp)
# target_link_libraries(SimpleTest ${libs})
#
# add_executable(VerticesATLAS main-vertices-atlas.cpp)
# target_link_libraries(VerticesATLAS ${plot-libs})
#
# add_executable(VerticesCMS main-vertices-cms.cpp)
# target_link_libraries(VerticesCMS ${plot-libs})
#
# add_executable(VerticesLHCb main-vertices-lhcb.cpp)
# target_link_libraries(VerticesLHCb ${plot-libs})
#
# add_executable(Conservation main-conservation.cpp)
# target_link_libraries(Conservation ${hepmc-libs})

# macro(test name file)
# add_executable(${name} ${file})
# target_link_libraries(${name} ${hepmc-libs})
# add_test(NAME ${name} COMMAND ${name})
# endmacro()

# test(detector test.cpp)
