/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-05-25 01:49:26
 */

#ifndef REVISELIST_H
#define REVISELIST_H

#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include "info.h"
#include "stringOpt.h"

using namespace std;
typedef pair<string, string> str2str;
struct Revision {
  vector<str2str> chglist;

  Revision(const string &);
  bool empty() const;
  int revise(string &);
  int revise(vector<string> &);
};

#endif
