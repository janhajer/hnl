cmake_minimum_required (VERSION 2.8)

project(Neutrino)

set (CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native -Ofast -Wall -Wextra -Wpedantic -Wfatal-errors -Wno-language-extension-token")

find_package(ROOT)
find_package(Boost COMPONENTS system filesystem REQUIRED)

include_directories("${ROOT_INCLUDE_DIRS}")
include_directories("${Boost_INCLUDE_DIR}")

link_directories(~/local/lib/)
add_library(libDelphes SHARED IMPORTED)

add_executable (CMS13 main-cms-13.cpp)
target_link_libraries(CMS13 Delphes ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

add_executable (CMS14 main-cms-14.cpp)
target_link_libraries(CMS14 Delphes ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

add_executable (LHCb main-lhcb.cpp)
target_link_libraries(LHCb Delphes ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

add_executable (Distance main-distance.cpp)
target_link_libraries(Distance Delphes ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

add_executable (Test main-test.cpp)
target_link_libraries(Test Delphes ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

add_executable (Test2 main-ray.cpp)
target_link_libraries(Test2 Delphes ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

add_executable (Wboson main-W.cpp)
target_link_libraries(Wboson Delphes ${ROOT_LIBRARIES} ${Boost_LIBRARIES})

add_executable (Leading main-lead.cpp)
target_link_libraries(Leading Delphes ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
