###
# Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences
# See the accompanying Manual for the contributors and the way to cite this work.
# Comments and suggestions welcome. Please contact
# Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
# 
# @Author: Dr. Guanghong Zuo
# @Date: 2020-12-07 09:06:29
# @Last Modified By: Dr. Guanghong Zuo
# @Last Modified Time: 2025-04-24 Thursday 18:00:30
###


## Stage for build cvtree
FROM alpine AS dev
LABEL Version=0.1 \
  MAINTAINER="Guanghong Zuo<ghzuo@ucas.ac.cn>"\
  description="Docker image for Collapse" 

## for develop environment
RUN apk --update add --no-cache g++ make cmake zlib-dev zlib-static
RUN apk --update add --no-cache nlohmann-json --repository=http://dl-cdn.alpinelinux.org/alpine/edge/community

## Build cltree
WORKDIR /root
COPY ../src /root/collapse
RUN mkdir collapse/build/ && cd collapse/build/ && cmake .. && make 

## Stage for run cvtree 
FROM alpine AS run
COPY --from=dev /root/collapse/build/bin/* /usr/local/bin/
RUN apk --update add --no-cache libstdc++

## for workplace
WORKDIR /root/data
ENTRYPOINT ["cltree"]
