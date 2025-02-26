cmake_minimum_required(VERSION 3.14)

include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(fcgr VERSION 1.0.0)
add_executable(fcgr main.cpp) 
set_target_properties(fcgr PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_compile_options(fcgr PRIVATE -Wall -Wextra)

# find_package(OpenMP REQUIRED)

# cnpy
message(STATUS "cnpy will be built from source")
FetchContent_Declare(cnpy
  GIT_REPOSITORY https://github.com/rogersce/cnpy.git
  GIT_TAG 4e8810b1a8637695171ed346ce68f6984e585ef4)
FetchContent_MakeAvailable(cnpy)

# KMC
message(STATUS "KMC3 will be built from source")
ExternalProject_Add(kmc
    GIT_REPOSITORY https://github.com/refresh-bio/KMC.git
    GIT_TAG d6c24dc88010508336d824bbffbf1ef3f49fddf8 # v3.2.4
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make kmc
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(kmc SOURCE_DIR)
SET(KMC_SOURCE_DIR ${SOURCE_DIR})
SET(KMC_INCLUDE_DIR ${KMC_SOURCE_DIR})
add_library(KMC STATIC IMPORTED)
set_target_properties(KMC PROPERTIES IMPORTED_LOCATION ${KMC_INCLUDE_DIR}/bin/libkmc_core.a)
add_dependencies(KMC kmc)

add_dependencies(fcgr cnpy KMC)

target_include_directories(fcgr
  PRIVATE ${KMC_INCLUDE_DIR}
  PRIVATE ${cnpy_SOURCE_DIR}
  PRIVATE ${PROJECT_SOURCE_DIR}
)

target_link_libraries(fcgr
  PUBLIC KMC
  PUBLIC z
  PUBLIC cnpy
)
