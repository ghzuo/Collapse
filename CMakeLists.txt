###
# Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai, China.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2020-12-07 09:06:29
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2020-12-07 09:37:48
###

cmake_minimum_required(VERSION 3.0)
#cmake_policy(SET CMP0012 NEW)
project (Collapse)
enable_language(CXX)

## For Compile
##########################################################
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
LINK_DIRECTORIES(${LIBRARY_OUTPUT_PATH})
SET(CMAKE_VERBOSE_MAKEFILE on)
ADD_DEFINITIONS("-O3 -std=c++11")
#set(CMAKE_FIND_LIBRARY_SUFFIXES .a .tbd .dylib .so)

## include subdirectories
ADD_SUBDIRECTORY(kit)
ADD_SUBDIRECTORY(collapse)



