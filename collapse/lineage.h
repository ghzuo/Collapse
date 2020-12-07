/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-07-21 11:48:01
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-07 14:28:28
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

#include "info.h"
#include "reviseList.h"

typedef unordered_map<string, string> TaxMap;
typedef pair<string, string> RankName;

// the lineage 
struct Lineage {
  bool def = true;
  string name;

  Lineage() = default;
  Lineage(const string& str) : name(str){};
};

struct LineageHandle {
  const string undefStr{"Unclassified"};
  regex misKingdom{"<D>([A-Za-z]+)<K>"+undefStr};
  regex prefix{"<[A-Za-z]>"};

  string taxfile;
  string revfile;
  map<string, char> rankmap{
      {"superkingdom", 'D'}, {"kingdom", 'K'},     {"phylum", 'P'},
      {"subphylum", 'p'},    {"superphylum", 'h'}, {"class", 'C'},
      {"subclass", 'c'},     {"superClass", 'l'},  {"order", 'O'},
      {"suborder", 'o'},     {"superorder", 'r'},  {"family", 'F'},
      {"subfamily", 'f'},    {"superfamly", 'a'},  {"genus", 'G'},
      {"subgenus", 'g'},     {"supergenus", 'e'},  {"species", 'S'},
      {"subspecies", 's'},   {"varietas", 'V'},    {"subvariety", 'v'},
      {"tribe", 'B'},        {"subtribe", 'b'},    {"section", 'n'},
      {"serotype", 'Y'},     {"isolate", 't'},     {"clade", 'd'}};
  vector<pair<string, char>> outrank;
  string outRankStr{"DKPCOFGS"};

  /// basic setting options
  LineageHandle() = default;
  LineageHandle(const string &);
  LineageHandle(const string &, const string &);
  LineageHandle(const string &, const string &, const string &, const string&);

  // set lineage handle
  void setRankByFile(const string &);
  void setOutRank();
  void setOutRankByTaxfile();
  size_t nOutRanks() const;

  // search entry
  void getLineage(vector<Lineage> &);
  void format(string &);
  void readTaxFile(TaxMap &);
  void checkRepeats(vector<Lineage> &);

  // convert linage style
  string lineageString(const RankName &) const;
  string lineageString(const vector<RankName> &) const;
};

size_t parseLineage(const string &, vector<string> &);
size_t separateLineage(const string &, vector<string> &);
string lastName(const string &);
string _lastName(const string &);
string commonLineage(const string &, const string &);
string commonLineage(const vector<string>&);
size_t nRanks(const string&);

#endif
