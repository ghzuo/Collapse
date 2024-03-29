/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-04-11 14:34:29
 */

#ifndef TAXARANK_H
#define TAXARANK_H

#include "kit.h"
#include <map>
#include <regex>
#include <string>
using namespace std;

typedef pair<string, string> RankName;

struct TaxaRank {
  // for keep one copy in program
  static TaxaRank *pRank;
  static TaxaRank *create();

  // mark option
  static pair<char, char> mark;
  template <class T> static string addMark(T sym) {
    string title;
    title += mark.first;
    title += sym;
    title += mark.second;
    return title;
  }

  // symbols for option
  string undefStr{"Unclassified"};
  string undefSym{mark.second + undefStr + mark.first};
  string strainMark{addMark('T')};
  regex misKingdom{addMark('D') + "([A-Za-z]+)" + addMark('K') + undefStr};
  string comKingdom{addMark('D') + "$1" + addMark('K') + "$1"};

  // taxon rank map and list
  map<string, char> rankmap{
      {"superkingdom", 'D'}, {"kingdom", 'K'},  {"subkingdom", 'k'},
      {"superphylum", 'Q'},  {"phylum", 'P'},   {"subphylum", 'p'},
      {"superclass", 'L'},   {"class", 'C'},    {"subclass", 'c'},
      {"superorder", 'W'},   {"order", 'O'},    {"suborder", 'o'},
      {"superfamly", 'M'},   {"family", 'F'},   {"subfamily", 'f'},
      {"genus", 'G'},        {"subgenus", 'g'}, {"species", 'S'},
      {"subspecies", 's'},   {"varietas", 'V'}, {"subvariety", 'v'},
      {"tribe", 'R'},        {"subtribe", 'r'}, {"section", 'E'},
      {"subsection", 'E'},   {"serotype", 'Y'}, {"isolate", 'I'},
      {"infraorder", 'i'},   {"biotype", 'B'},  {"genotype", 'N'},
      {"domain", 'D'}};
  vector<pair<string, char>> outrank{
      {"Domain", 'D'}, {"Kingdom", 'K'}, {"Phylum", 'P'}, {"Class", 'C'},
      {"Order", 'O'},  {"Family", 'F'},  {"Genus", 'G'},  {"Species", 'S'}};
  bool fixed{false};

  // set rank handle
  void initial(const string &, const string &);
  string setRankByFile(const string &);
  void setOutRank(const string &, bool fix = false);
  void setOutRank(const vector<pair<string, char>> &, bool fix = false);

  // output rank in formats
  void outRanksJson(ostream &);
  void outRanksCSV(ostream &);

  // for check
  size_t nOutRanks() const;
  bool wellDefined(const string &);
  void format(string &);
  string lineage(const string &, const string &);

  // convert linage style
  char getSymbol(const string &) const;
  string rankString(const RankName &) const;
  string rankString(const vector<RankName> &) const;

  //output the rank index
  int rankindex(const string&);

private:
  TaxaRank() = default;
};

// comman options
size_t parseLineage(const string &, vector<string> &);
size_t separateLineage(const string &, vector<string> &);
string lastName(const string &);
string lastNameNoRank(const string &);
string commonLineage(const string &, const string &);
string commonLineage(const vector<string> &);
string delStrain(const string &);
size_t nRanks(const string &);

#endif // !TAXARANK_H