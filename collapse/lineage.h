/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-07-21 11:48:01
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-28 12:10:53
 */

#ifndef LINEAGE_H
#define LINEAGE_H

#include <cctype>
#include <fstream>
#include <map>
#include <regex>
#include <string>
#include <unordered_map>
#include <vector>

#include "kit.h"
#include "reviseList.h"
#include "taxarank.h"
#include "taxadb.h"

typedef unordered_map<string, string> TaxMap;

// the lineage
struct Lineage {
  bool def = false;
  string name;

  Lineage() = default;
  Lineage(const string &str) : name(str){};
  friend ostream& operator<<(ostream&, const Lineage&);
};

struct LineageHandle { 
  string tdbpath; 
  string taxfile;
  string revfile;
  TaxaRank* rank;
  
  /// basic setting options
  LineageHandle();
  LineageHandle(const string &, const string &, const string &);

  // search entry
  void getLineage(vector<Lineage> &);
  void readTaxFile(TaxMap &);
  void checkRepeats(vector<Lineage> &);
};

#endif
