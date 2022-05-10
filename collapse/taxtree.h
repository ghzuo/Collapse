/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-05-10 10:51:01
 */

#ifndef TREE_H
#define TREE_H

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "kit.h"
#include "taxarank.h"
using namespace std;

typedef pair<string, string> str2str;

struct Node {
  typedef vector<Node *> Children;

  string name;
  size_t id;
  double length;
  double bootstrap;
  double depth;
  Node *parent;
  Children children;

  size_t taxSize, taxLevel, nleaf, nxleaf,
      nUpLeaf; // nxleaf is the number of unclassfied leafs
  bool unclassified, uploaded, otu;

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
  void swap(Node *);

  void checkUnclassified();
  void checkUploaded();
  void getDefineLeafs(vector<Node *> &);
  void getUndefineLeafs(vector<Node *> &);
  void getUndefineNames(vector<string> &);
  void chkLeafsName(vector<string> &, vector<string> &);

  void setOneLeaf(const string &, bool);
  void _setOneBranch();
  void setAllBranches();
  void setBranchLineage();
  string _getBranchName(const vector<Node *> &);
  int nClade();

  void _getPrediction(string &);
  void outPrediction(ostream &);
  void outPrediction(const string &);

  Node *resetroot(const string &);
  Node *resetroot(Node *);

  void updateRootedTree();
  Node *rootingDirect();
  Node *rootingByOutgrp(const string &);
  Node *rootingByTaxa();
  bool _findOutgrpCandidates(vector<Node *> &);
  Node *_forceRooting(Node *);

  void balanceTree(const string &, size_t);
  void findRootCandidates(vector<Node *> &, size_t);
  void _findRootCandidates(vector<Node *> &, size_t);
  void _rearrangeOutgroup(Node *);
  Node *_sibling();

  void _madTree(const vector<Node *> &);
  pair<double, double> getMAD();
  void _getOTUad(vector<double>&, double&);

  void _pmrTree(const vector<Node *> &);
  double getPMR();
  void _getOTUdepth(vector<double> &);

  void _mdmpTree(const vector<Node *> &);
  void _setLengthByMidpoint();
  void _getDepth();

  void _mpTree();
  void _getMaxPath(pair<double, vector<Node*>>&, pair<double, vector<Node*>>&);

  void _outnwk(ostream &);
  static function<string(Node *)> nwkname;
  void outnwk(ostream &);
  void outnwk(const string &);
  void outitol(const string &);
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
  void updateId(int &);
  void reinitTree();
  void chgLeafName(const str2str &);
};

#endif
