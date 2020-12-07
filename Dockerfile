###
# Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai, China.
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@fudan.edu.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2020-12-07 09:06:29
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2020-12-07 09:56:59
###

## Stage for build cvtree
FROM alpine AS dev
LABEL Version=0.1 \
  MAINTAINER="Guanghong Zuo<ghzuo@fudan.edu.cn>"\
  description="Docker image for Collapse" 

## for develop environment
RUN apk --update add --no-cache g++ make cmake zlib-dev

## Build cvtree
WORKDIR /root
COPY ./collapse /root/collapse/collapse
COPY ./kit /root/collapse/kit
COPY ./CMakeLists.txt /root/collapse/
RUN mkdir collapse/build/ && cd collapse/build/ && cmake .. && make 

## Stage for run cvtree 
FROM alpine AS run
COPY --from=dev /root/collapse/build/bin/* /usr/local/bin/
RUN apk --update add --no-cache libstdc++

## for workplace
WORKDIR /root/data
