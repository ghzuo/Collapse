/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-12-05 15:05:30
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-06 15:57:20
 */

#ifndef DMP2LINEAGE_H
#define DMP2LINEAGE_H
#include <algorithm>
#include <regex>
#include <string>
#include <vector>

#include "info.h"
#include "lineage.h"
#include "stringOpt.h"

using namespace std;

struct TaxonNode {
  size_t parent = 0;
  size_t taxid = 0;
  string name;
  string level;
  string lineage;

  friend ostream &operator<<(ostream &, const TaxonNode &);
};

struct TaxaDB {
  vector<TaxonNode> taxlist;
  LineageHandle hdl;

  // the regex for name
  regex nonWord{"[^A-z0-9]+"};
  regex regTaxID{"taxid([0-9]+)"};
  regex regSpeciesName{"^([A-Z][a-z]+_[a-z]+)"};
  regex regGenusName{"^([A-Z][a-z]+)"};

  TaxaDB() = default;

  // build from dump files
  TaxaDB(const string&, const string&);
  void _readNodeDump(const string &);
  void _readNameDump(const string &);
  
  // set lineage map
  void resetRankMap(const string&);

  //get lineage
  string search(const string&);
  TaxonNode getLineage(size_t);
  TaxonNode getLineage(const string&);
  void exportLineage(ostream& os);
  string _getlng(size_t);
};

#endif // !DMP2LINEAGE_H
