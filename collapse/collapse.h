/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-09-01 15:53:05
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2021-01-18 14:17:41
 */

#ifndef COLLAPSE_H
#define COLLAPSE_H

#include <fstream>
#include <sstream>

#include "kit.h"
#include "taxsys.h"
#include "taxtree.h"
using namespace std;

// read arguments
struct RunArgs {
  string program;
  string infile;
  string taxadb, taxfile, taxrev;
  string abfile, abtype;
  string outPref;
  string outgrp;
  bool forWeb, forApp, predict; 
  // two hidden options for output for server and app

  RunArgs(int, char **);
  void usage();
};

void collapse(int, char**);

void output(const vector<Lineage>&, Taxa &, Node *, RunArgs &);
void out4serv(const vector<Lineage>&, Taxa&, Node *, RunArgs &);
void out4app(const vector<Lineage>&, Taxa&, Node *, const string&);
void outTaxaJson(Taxa&, Node *, ostream&);
void outTreeJson(Taxa&, Node *, ostream&);
void outLngsJson(const vector<Lineage>&, Taxa&, ostream&);
#endif
