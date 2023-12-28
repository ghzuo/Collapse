###
# Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2021-01-22 15:47:25
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2022-03-16 16:41:51
###

cmake_minimum_required(VERSION 3.6)
#cmake_policy(SET CMP0012 NEW)
project (CLTree)
enable_language(CXX)

## For Compile
##########################################################
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
LINK_DIRECTORIES(${LIBRARY_OUTPUT_PATH})
# SET(CMAKE_CXX_FLAGS "-Wall")
SET(CMAKE_VERBOSE_MAKEFILE on)
ADD_DEFINITIONS("-O3 -std=c++11")
#set(CMAKE_FIND_LIBRARY_SUFFIXES .a .tbd .dylib .so)

### use conda library
if (NOT $ENV{CONDA_PREFIX} STREQUAL "")
  message(STATUS "CURRENT CONDA PREFIX: $ENV{CONDA_PREFIX}")
  include_directories(BEFORE $ENV{CONDA_PREFIX}/include)
  #link_directories(BEFORE $ENV{CONDA_PREFIX}/lib)
endif ()

### Check required LIABRARIES
if(NOT EMSCRIPTEN)
  SET(ZLIB_USE_STATIC_LIBS ON)
  find_package(ZLIB)
  if(ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    link_libraries(${ZLIB_LIBRARIES})
  else()
    message(FATAL_ERROR "Zlib not found")
  endif()
endif()

### find the nlohmann json
find_package(nlohmann_json 3.6 REQUIRED)

## include subdirectories
ADD_SUBDIRECTORY(kit)
ADD_SUBDIRECTORY(collapse)



