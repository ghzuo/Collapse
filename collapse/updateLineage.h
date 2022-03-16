/*
 * Copyright (c) 2021  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2021-06-21 14:40:41
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2021-06-21 15:08:11
 */

#ifndef UPDATELINEAGE_H
#define UPDATELINEAGE_H

#include <fstream>
#include <sstream>

#include "kit.h"
#include "lineage.h"
using namespace std;

// read arguments
struct UpLngArgs {
  string program;
  string taxadb, taxfile, taxrev;
  string abfile, abtype;
  vector<string> nmlist;
  string outfile;
  string format;
  // two hidden options for output for server and app

  UpLngArgs(int, char **);
  void usage();
};

void updateLineage(int, char**);
#endif
