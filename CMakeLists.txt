cmake_minimum_required(VERSION 3.1)
project(OLSQ
VERSION 0.1.0
LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_BUILD_TYPE release)
# set(CMAKE_BUILD_TYPE Debug)

############################################################
# Bitwuzla
############################################################

find_library(BITWUZLA_LIBRARY bitwuzla REQUIRED)
find_library(BITWUZLABV_LIBRARY bitwuzlabv REQUIRED)
find_library(BITWUZLABB_LIBRARY bitwuzlabb REQUIRED)
find_library(BITWUZLALS_LIBRARY bitwuzlals REQUIRED)
get_filename_component(BITWUZLA_ROOT ${BITWUZLA_LIBRARY} DIRECTORY)
get_filename_component(BITWUZLA_ROOT ${BITWUZLA_ROOT} DIRECTORY)
get_filename_component(BITWUZLA_ROOT ${BITWUZLA_ROOT} DIRECTORY)
set(BITWUZLA_INCLUDE_DIR "${BITWUZLA_ROOT}/include")

############################################################
# gmp
############################################################
find_library(GMP_LIBRARY gmp REQUIRED)
get_filename_component(GMP_ROOT ${GMP_LIBRARY} DIRECTORY)
set(GMP_INCLUDE_DIR "${GMP_ROOT}/include")
# message(STATUS "Find gmp: ${GMP_ROOT}")

############################################################
# OLSQ library target
############################################################


add_library(olsq
  src/cir/circuit.cpp
  src/device/device.cpp
  src/misc/timeUsage.cpp
  src/clusterer/clusterer2.cpp
  src/clusterer/clusterer3.cpp
  src/placer/initialMapper.cpp
  src/molsq/molsq.cpp
  src/olsq2/olsq2.cpp
  src/stbolsq2/stbolsq2.cpp
  src/api/apiPy.cpp
  src/placer/placer.cpp
  src/placer/placer_dev.cpp
  src/router/aRouter.cpp
  src/router/aRouter_dev.cpp 
  src/writer/writer.cpp
  src/rolsq2/rolsq2.cpp
)

set(SRC
  src/cir/circuit.cpp
  src/device/device.cpp
  src/misc/timeUsage.cpp
  src/clusterer/clusterer2.cpp
  src/clusterer/clusterer3.cpp
  src/placer/initialMapper.cpp
  src/molsq/molsq.cpp
  src/olsq2/olsq2.cpp
  src/stbolsq2/stbolsq2.cpp
  src/api/apiPy.cpp
  src/placer/placer.cpp
  src/placer/placer_dev.cpp
  src/router/aRouter.cpp
  src/router/aRouter_dev.cpp 
  src/writer/writer.cpp
  src/rolsq2/rolsq2.cpp
)

set(API
   src/api/apiPy.cpp
 )

target_include_directories(olsq PUBLIC src include
  "$<BUILD_INTERFACE:${BITWUZLA_INCLUDE_DIR}>"
  "$<BUILD_INTERFACE:${GMP_INCLUDE_DIR}>"
)

set_property(SOURCE src/router/aRouter_dev.cpp
  PROPERTY COMPILE_FLAGS "-Wall -Werror -Wpedantic -Wno-error=unused-variable -Wno-error=format"
)

find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})


############################################################
# Install Boost                
############################################################
SET(BASEPATH ${PROJECT_SOURCE_DIR})
#SET(BASEPATH .)
SET(SRCPATH ${BASEPATH}/src)
SET(BINPATH ${BASEPATH}/bin)
SET(INCPATH ${BASEPATH}/include)
SET(LIBPATH ${BASEPATH}/lib)

IF (NOT EXISTS "${INCPATH}")
  FILE(MAKE_DIRECTORY ${INCPATH})
ENDIF()

find_package(Boost)
IF (NOT EXISTS "${LIBPATH}")
  FILE(MAKE_DIRECTORY ${LIBPATH})
ENDIF()

INCLUDE_DIRECTORIES(${SRCPATH} ${INCPATH})

# Install Boost
SET(BOOST_VERSION 1_82_0)
FIND_PACKAGE(Boost 1.82)
IF(Boost_FOUND)
  SET(Boost_USE_STATIC_LIBS   ON)
  SET(Boost_INCLUDE_DIR ${INCPATH}/boost_${BOOST_VERSION})
  SET(Boost_LIBRARY_DIR_RELEASE ${INCPATH}/boost_${BOOST_VERSION})
  SET(Boost_USE_STATIC_LIBS   ON)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ELSE()
  IF (NOT EXISTS "${INCPATH}/boost_${BOOST_VERSION}")
    execute_process (COMMAND bash -c "wget https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.gz -P ${INCPATH}")
    #execute_process (COMMAND bash -c "wget https://dl.bintray.com/boostorg/release/1.72.0/source/boost_1_72_0.tar.gz -P ${INCPATH}")
    execute_process (COMMAND bash -c "tar zxvf ${INCPATH}/boost_${BOOST_VERSION}.tar.gz -C ${INCPATH}")
    execute_process (COMMAND bash -c "rm ${INCPATH}/boost_${BOOST_VERSION}.tar.gz")
  ENDIF()

  SET(BOOST_ROOT ${INCPATH}/boost_${BOOST_VERSION})
  SET(Boost_INCLUDE_DIR ${INCPATH}/boost_${BOOST_VERSION})
  SET(Boost_LIBRARY_DIR_RELEASE ${INCPATH}/boost_${BOOST_VERSION})
  SET(Boost_USE_STATIC_LIBS   ON)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF()

############################################################
# library                         
############################################################
set(LIBRARY_OUTPUT_PATH ${LIBPATH})
# LIB
file(GLOB LIB ${LIBPATH}/*.a)


set(PYBIND11_ROOT_DIR "include/pybind11")
message(STATUS "Find Pybind11: ${PYBIND11_ROOT_DIR}")

include_directories("${PYBIND11_ROOT_DIR}/include")
add_subdirectory("include/pybind11")
set(PYBIND11_CPP_STANDARD -std=c++14)

#find_package(Bitwuzla REQUIRED)
include_directories(include/pblib/include/pblib)


target_link_libraries(olsq PUBLIC ${PROJECT_SOURCE_DIR}/include/pblib/lib/libpb.a ${GMP_LIBRARY} ${BITWUZLA_LIBRARY} ${BITWUZLABB_LIBRARY} ${BITWUZLABV_LIBRARY} ${BITWUZLALS_LIBRARY} ${CMAKE_DL_LIBS})

add_executable(main app/main.cpp)
target_link_libraries(main olsq ${CMAKE_DL_LIBS})



pybind11_add_module(olsqPy ${SRC} ${API})
target_include_directories(olsqPy PUBLIC 
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>"
  "$<BUILD_INTERFACE:${BITWUZLA_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(olsqPy PUBLIC ${PROJECT_SOURCE_DIR}/include/pblib/lib/libpb.a ${GMP_LIBRARY} ${BITWUZLA_LIBRARY} ${BITWUZLABB_LIBRARY} ${BITWUZLABV_LIBRARY} ${BITWUZLALS_LIBRARY} ${CMAKE_DL_LIBS})

# find_package(Python3 REQUIRED COMPONENTS Development)
# find_package(PythonLibs 3)

