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
set(HEADERS
        src/reco_acceptance.h
        src/sim_acceptance.h
        )

set(SOURCES
        src/reco_acceptance.cc
        src/sim_acceptance.cc
        )
add_definitions(${ROOT_CXX_FLAGS})

add_executable(acceptance src/acceptance.cc ${SOURCES} ${HEADERS})
target_link_libraries(acceptance ${ROOT_LIBRARIES} ROOT::EG AnalysisTreeBase AnalysisTreeInfra)

add_executable(efficiency src/efficiency.cc)
target_link_libraries(efficiency ${ROOT_LIBRARIES})

add_executable(fit_layers src/fit_layers.cc)
target_link_libraries(fit_layers ${ROOT_LIBRARIES})

add_executable(fit_occupancy src/fit_occupancy.cc)
target_link_libraries(fit_occupancy ${ROOT_LIBRARIES})

add_executable(equal_bins src/cut_equal_bins.cc)
target_link_libraries(equal_bins ${ROOT_LIBRARIES})

add_executable(theta_p src/theta_p.cc)
target_link_libraries(theta_p ${ROOT_LIBRARIES})

add_executable(occupancy src/calculate_occupancy.cc)
target_link_libraries(occupancy ${ROOT_LIBRARIES})