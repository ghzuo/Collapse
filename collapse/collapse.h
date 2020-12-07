/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-09-01 15:53:05
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-11-27 21:40:49
 */

#ifndef COLLAPSE_H
#define COLLAPSE_H

#include <fstream>
#include <sstream>

#include "info.h"
#include "reviseList.h"
#include "taxsys.h"
#include "taxtree.h"
using namespace std;

// read arguments
struct Args {
  string program;
  string infile;
  string taxfile, taxrev;
  string abfile, abtype;
  string outPref, treeSuff;
  string outgrp;
  bool forWeb, predict;

  Args(int, char **);
  void usage();
};

void output(Taxa &, Node *, Args &);
#endif
