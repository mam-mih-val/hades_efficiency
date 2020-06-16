cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(efficiency)
set(CMAKE_CXX_STANDARD 17)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH $ENV{AnalysisTree_HOME})
list(APPEND CMAKE_PREFIX_PATH ${AnalysisTree_HOME})

find_package(ROOT REQUIRED COMPONENTS RIO Tree)
find_package(AnalysisTree REQUIRED)

set(CMAKE_CXX_FLAGS "-Wall")

include(${ROOT_USE_FILE})

include_directories(
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
        ${ROOT_INCLUDE_DIRS}
        ${AnalysisTree_INCLUDE_DIR}
)

link_directories(
        $ENV{HADDIR}/lib
        ${CMAKE_SOURCE_DIR}
        ${AnalysisTree_LIBRARY_DIR}
)

set(SOURCES
        src/momentum_reconstruction_eff.cc
        )
add_definitions(${ROOT_CXX_FLAGS})

add_executable(efficiency src/main.cc ${HEADERS} ${SOURCES})
target_link_libraries(efficiency ${ROOT_LIBRARIES} AnalysisTreeBase)