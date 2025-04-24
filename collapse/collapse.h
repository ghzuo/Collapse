/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2025-03-16 Sunday 12:07:17
 */

#ifndef COLLAPSE_H
#define COLLAPSE_H

#include <fstream>
#include <sstream>
#include <set>

#include "kit.h"
#include "taxsys.h"
#include "taxtree.h"
using namespace std;

// read arguments
struct RunArgs {
  string program;
  string infile;
  string taxadb, taxfile, taxrev;
  string rankfile, outrank;
  string outPref;
  string lngfile;
  string outgrp;
  string clevel;
  string rootMeth;
  string otuLevel;
  bool forWeb, forApp, predict;
  bool itol, byBranch;
  // two hidden options for output for server and app

  RunArgs(int, char **);
  void usage();
};

void collapse(int, char **);

void output(const LngData &, Taxa &, Node *, RunArgs &);
void out4serv(const LngData &, Taxa &, Node *, RunArgs &);
void out4app(const LngData &, Taxa &, Node *, const string &);
void outTaxaJson(Taxa &, Node *, ostream &);
void outTreeJson(Taxa &, Node *, ostream &);
void outLngsJson(const LngData &, Taxa &, ostream &);

void outItolNodes(const vector<Node *>&, const string &);
void outItolLabel(const vector<Node *>&, const string &);
void outItolPopup(const Taxa &, const vector<Node *>&, const string &);
void outItolSymbol(const Taxa &, const vector<Node *>&, const string &);
void outItolStrap(const vector<Node *>&, const map<string, TaxonState>&, const string &);
void outItolCollapse(const vector<Node*>&, const map<string, TaxonState>&, const string&);
void itolHeader(ostream&, const string&, const string&);
void getDivision(const Taxa&, const string&, map<string, TaxonState>&);
void getTopDivision(const Taxa &, map<string, TaxonState>&);
void getColorMap(map<string,string>&);
string itolPopusStr(Node*, const Taxa&);

#endif
