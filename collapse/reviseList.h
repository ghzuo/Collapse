/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-07-21 11:48:01
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-05 09:06:11
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
  void revise(string &);
  void revise(vector<string> &);
};

#endif
