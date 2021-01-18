/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-01-03 21:03:33
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-21 18:10:51
 */

#ifndef TAXSYS_H
#define TAXSYS_H

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <stdio.h>
#include <string>
#include <tuple>
#include <vector>

#include "kit.h"
#include "lineage.h"
#include "taxtree.h"

using namespace std;

typedef vector<pair<string, char>> Abbr;

/********************************************************************************
 * @brief to record the state of one taxon level, e.g. number of monophyly,
 * entropy
 *
 ********************************************************************************/
struct TaxLevelState {
  int nSolo, nMono, nPoly;
  double sTax, sTree;
  TaxLevelState() : nSolo(0), nMono(0), nPoly(0), sTax(0), sTree(0){};
};

/********************************************************************************
 * @brief to record state of a taxon
 *
 ********************************************************************************/
struct TaxonState {
  bool monophy;
  size_t nStrain;
  vector<size_t> distract;
  TaxonState() : monophy(false), nStrain(1){};
};

/********************************************************************************
 * @brief the taxonomy system for statistics and annotate the tree
 *
 ********************************************************************************/
struct TaxSys {
  const static string rootTaxon;

  size_t nStrain;
  map<string, TaxonState> state;

  TaxSys() = default;
  TaxSys(const vector<string> &);
  void initial(const vector<string> &);

  void annotateBranch(size_t, int, string &, size_t &);
  void annotateLeaf(int, string &, size_t &);
};

/********************************************************************************
 * @brief set of two taxsyses of defined lineage and undefined lineage (with
 * Unclassified items)
 ********************************************************************************/
struct Taxa {
  TaxSys def, undef;
  TaxaRank *rank;
  Taxa(const vector<Lineage> &);

  void annotate(Node *);

  void outTax(ostream &);
  void outTax(const string &);
  void outUnclass(vector<string> &, ostream &);
  void outUnclass(vector<string> &, const string &);
  void outStatitics(ostream &);
  void outStatitics(const string &);
  void outEntropy(ostream &);
  void outEntropy(const string &);
  void outJsonEntropy(ostream &);
  void outJsonTax(ostream &);
  void outJsonUnclass(vector<string> &, ostream &);
  void outJson(const string &);
};

#endif
