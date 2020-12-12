/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-12-05 15:05:30
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-12 09:59:23
 */

#ifndef TAXADB_H
#define TAXADB_H
#include <algorithm>
#include <regex>
#include <string>
#include <vector>

#include "kit.h"
#include "taxarank.h"

using namespace std;
/********************************************************************************
 * @brief class for a taxon node
 *
 ********************************************************************************/
struct TaxonNode {
  size_t parent = 0;
  size_t taxid = 0;
  string level;
  string name;
  vector<string> nicks;
  string lineage;

  friend ostream &operator<<(ostream &, const TaxonNode &);
};

/********************************************************************************
 * @brief class for the taxonomy system
 *
 ********************************************************************************/
struct TaxaDB {
  vector<TaxonNode> taxlist;
  map<string, size_t> nmlist;
  TaxaRank* rank;

  // the regex for name
  regex nonWord{"[^A-Za-z0-9]+"};
  regex regTaxID{"taxid([0-9]+)"};


  // build from dump files
  TaxaDB(const string &);
  string goodname(const string&);
  string idname(const string&);
  void tgz4taxdb(const string&);
  void _readNodeDump(istream &);
  void _readNameDump(istream &);

  // read/write gzip database file
  void readTable(const string &);
  void writeTable(const string &);

  // search lineage
  string search(size_t);
  string search(const string &);
  TaxonNode getLineage(size_t);
  size_t _name2id(const string&);
  size_t _candidateQuery(const string&, vector<string>&);
  string _getlng(size_t);
  void _initNameMap();

  // query full lineage node
  void query(size_t, vector<TaxonNode>&);
  void query(const string&, vector<TaxonNode>&);
  bool getFullLineage(size_t, vector<TaxonNode>&);
  void _reQuery(size_t, vector<TaxonNode>&);
};

#endif // !TAXADB_H
