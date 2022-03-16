/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-16 12:29:26
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
#include "taxadb.h"
#include "taxarank.h"

// the lineage
struct Lineage {
  bool def = false;
  string name;

  Lineage() = default;
  Lineage(const string &str) : name(str){};
  friend ostream &operator<<(ostream &, const Lineage &);
};

typedef unordered_map<string, string> TaxMap;
struct LngData {
  string tdbpath;
  string revfile;
	vector<string> tflist;
  vector<Lineage> data;
  TaxaRank *rank;

  /// basic setting options
  LngData();
  LngData(const string &, const string &, const string &);
  void initData(const vector<string> &);

  // search entry
  void getLineage(vector<string> &);
  size_t getLngFromFile();
	void getLngFromDB();
  void readListFile(istream &, TaxMap &);
  void readJsonFile(istream &, TaxMap &);
  void readCsvFile(istream &, TaxMap &);
  void checkRepeats();

  // output data
  Lineage &operator[](size_t);
  void getStat(vector<pair<int, int>>&) const;
	size_t size() const;

	void output(const string&, const string&) const;
  void outjson(ostream &) const;
	void outcsv(ostream &) const;
	void outlist(ostream &) const;
};

#endif
