/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-06-01 16:20:31
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-06-01 16:37:40
 */

#ifndef ROOTING_H
#define ROOTING_H

#include <fstream>
#include <sstream>
#include <set>

#include "kit.h"
#include "taxtree.h"
using namespace std;

// read arguments
struct RootingArgs {
  string program;
  string infile;
  string outfile;
  string rootMeth;
  string outgrp;

  RootingArgs(int, char **);
  void usage();
};

void rooting(int, char **);

#endif