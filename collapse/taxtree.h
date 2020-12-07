/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-03-17 15:39:23
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-07 16:46:50
 */

#ifndef TREE_H
#define TREE_H

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "info.h"
#include "lineage.h"
#include "stringOpt.h"
using namespace std;

typedef pair<string, string> str2str;

struct Node {
  typedef vector<Node *> Children;

  string name;
  size_t id;
  double length;
  double bootstrap;
  Node *parent;
  Children children;

  size_t taxSize, taxLevel, nleaf, nxleaf,
      nUpLeaf; // nxleaf is the number of unclassfied leafs
  bool unclassified, uploaded;

  Node();
  Node(size_t);
  Node(size_t, const string &);
  Node(size_t, const vector<Node *> &);

  void clear();
  void addChild(Node *);
  void deleteChild(Node *);
  void getDescendants(vector<Node *> &);
  void getLeafs(vector<Node *> &);
  void getBranches(vector<Node *> &);
  void getAllNodes(vector<Node *> &);
  bool isLeaf();

  void checkUnclassified();
  void checkUploaded();
  void getDefineLeafs(vector<Node *> &);
  void getUndefineLeafs(vector<Node *> &);
  void getUndefineNames(vector<string> &);
  void chkLeafsName(vector<string> &, vector<string> &);

  void setOneLeaf(const Lineage&);
  void _setOneBranch();
  void setAllBranches();
  void setBranchLineage();
  string _getBranchName(const vector<Node *> &);
  int nClade();

  string getStrainName();
  void _getPrediction(string &);
  void outPrediction(ostream &);
  void outPrediction(const string &);

  Node *resetroot(const string &);
  Node *resetroot(Node *);
  void updateRootedTree();
  Node *rootingDirect();
  Node *rootingByOutgrp(const string &);
  Node *rootingByTaxa();
  bool _getOutgrpCandidates(vector<Node *> &);
  Node *_forceRooting(Node *);

  void _outnwk(ostream &);
  void outnwk(ostream &);
  void outnwk(const string &);
  void _innwk(istream &);
  void _nwkItem(const string &);
  void innwk(istream &);
  void innwk(const string &);

  void _injson(istream &);
  void _getStr(istream &, string &);
  void _getKeyValue(string &, string &);
  void injson(istream &);
  void injson(const string &);
  void outjson(ostream &);
  void outjson(const string &);
  void outjsonAbbr(ostream &);

  void renewId(const unordered_map<string, size_t> &);
  void reinitTree();
  void chgLeafName(const str2str &);
};

#endif
