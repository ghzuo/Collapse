/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-05-07 18:19:19
 */

#include "taxtree.h"
const size_t N_FORKS(2);

/*******************************************************************/
/********* Member Functions For Node Class *************************/
/*******************************************************************/
Node::Node()
    : name(""), id(0), length(NAN), depth(NAN), parent(NULL), taxSize(0),
      taxLevel(0), nleaf(0), nxleaf(0), unclassified(false), uploaded(false),
      otu(false) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n)
    : name(""), id(n), length(NAN), depth(NAN), parent(NULL), taxSize(0),
      taxLevel(0), nleaf(0), nxleaf(0), unclassified(false), uploaded(false),
      otu(false) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n, const string &str)
    : name(str), id(n), length(NAN), depth(NAN), parent(NULL), taxSize(0),
      taxLevel(0), nleaf(0), nxleaf(0), unclassified(false), uploaded(false),
      otu(false) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n, const vector<Node *> &vn)
    : name(""), id(n), length(NAN), depth(NAN), parent(NULL), children(vn),
      taxSize(0), nleaf(0), nxleaf(0), unclassified(false), uploaded(false),
      otu(false) {
  children.reserve(2);
};

bool Node::isLeaf() { return children.empty(); };

void Node::addChild(Node *nd) {
  children.emplace_back(nd);
  (*nd).parent = this;
};

void Node::deleteChild(Node *nd) {
  Children::iterator iter = find(children.begin(), children.end(), nd);
  children.erase(iter);
};

void Node::clear() {
  vector<Node *> nodes;
  getDescendants(nodes);
  vector<Node *>::iterator iter = nodes.begin();
  vector<Node *>::iterator iterEnd = nodes.end();
  for (; iter != iterEnd; ++iter)
    delete *iter;
  children.clear();
}

void Node::getDescendants(vector<Node *> &nds) {
  if (!isLeaf()) {
    for (const auto &nd : children) {
      nds.emplace_back(nd);
      (*nd).getDescendants(nds);
    }
  }
};

void Node::getAllNodes(vector<Node *> &nds) {
  nds.emplace_back(this);
  getDescendants(nds);
};

void Node::getLeafs(vector<Node *> &nodes) {
  if (isLeaf()) {
    nodes.emplace_back(this);
  } else {
    for (const auto &nd : children)
      nd->getLeafs(nodes);
  }
};

void Node::getBranches(vector<Node *> &nodes) {
  if (!isLeaf()) {
    nodes.emplace_back(this);
    for (const auto &nd : children)
      nd->getBranches(nodes);
  }
};

void Node::swap(Node *ndx) {
  // swap content of two node
  Node tmp = *this;
  *this = *ndx;
  *ndx = tmp;

  // update the parent of children
  for (auto nd : children)
    nd->parent = this;
  for (auto nd : ndx->children)
    nd->parent = ndx;
}

/********************************************************************************
 * @brief tools for set root of the unroot tree by simple method
 *
 * @param str
 * @return Node*
 ********************************************************************************/
void Node::updateRootedTree() {
  theInfo("This tree is a rooted treee, keep as it is");

  // set branch lineage: type, fullname,  taxLevel
  setAllBranches();
}

Node *Node::rootingDirect() {
  Node *theTree = _forceRooting(this);
  theTree->setAllBranches();
  return theTree;
};

Node *Node::rootingByOutgrp(const string &str) {
  theInfo("Rooting the tree by leaf: " + str);
  Node *theTree = resetroot(str);
  if (theTree == NULL) {
    theInfo("Cannot find leaf named as " + str);
    return theTree;
  }

  // rooting the tree
  theTree = _forceRooting(theTree);

  // set branch lineage: type, fullname,  taxLevel
  theTree->setAllBranches();
  return theTree;
};

/********************************************************************************
 * @brief rooting tree by taxa
 *
 * @return Node*
 ********************************************************************************/
Node *Node::rootingByTaxa() {

  // preset the branch attraction
  Node *theTree = this;
  theTree->setAllBranches();

  // find the candidate of outgroup
  vector<Node *> clades;
  for (auto &nd : theTree->children) {
    nd->_findOutgrpCandidates(clades);
  }

  // root the tree by the highest rank clade and length
  if (!clades.empty()) {

    // preselect the outgroup as the smallest branch
    Node *outgrp =
        *(min_element(clades.begin(), clades.end(),
                      [](Node *a, Node *b) { return a->nleaf < b->nleaf; }));
    // record the change branch
    Node *chgBranch = outgrp->parent->parent;
    // rearrange the tree (still unroot tree)
    theTree = resetroot(outgrp);
    // renew the update branch
    chgBranch->setAllBranches();

    // find a good outgroup (the last item of children)
    // by higest rank of common lineage or longest branch length
    size_t minComLev = nRanks(commonLineage(chgBranch->name, outgrp->name));
    auto outIter = theTree->children.rbegin();
    for (auto iter = theTree->children.rbegin() + 1;
         iter != theTree->children.rend(); ++iter) {
      if (*iter != chgBranch) {
        size_t comLev = nRanks(commonLineage(chgBranch->name, (*iter)->name));
        if (comLev < minComLev) {
          outIter = iter;
          minComLev = comLev;
        }
      }
    }

    // swap new outgroup to the back
    if (outgrp != *outIter) {
      auto tmp = *outIter;
      *outIter = outgrp;
      theTree->children.back() = tmp;
    }
  }

  // rooting the unroot tree to a rooted tree
  theTree = _forceRooting(theTree);

  // renew the added node and the root
  if (!theTree->children.front()->isLeaf())
    theTree->children.front()->_setOneBranch();
  if (!theTree->children.back()->isLeaf())
    theTree->children.back()->_setOneBranch();
  theTree->_setOneBranch();

  return theTree;
};

bool Node::_findOutgrpCandidates(vector<Node *> &clades) {
  if (unclassified) {
    return false;
  } else if (isLeaf() || nClade() > 0) {
    return true;
  } else {
    // the candidates are the nodes whose all siblings are clades
    bool isCandidate(true);
    for (auto nd : children) {
      isCandidate = nd->_findOutgrpCandidates(clades) && isCandidate;
    }
    if (isCandidate) {
      clades.insert(clades.end(), children.begin(), children.end());
    }
    return false;
  }
};

Node *Node::_forceRooting(Node *root) {

  // testing whether a rooted tree
  if (root->children.size() == 2) {
    theInfo("This is a rooted tree");
    return root;
  }

  // add a super root
  Node *theRoot = new Node();

  // add the last child of node as the outgroup of theRoot
  // add root to the super root (theRoot)
  Node *outgrp = root->children.back();
  theRoot->addChild(outgrp); // outgrp at front
  theRoot->addChild(root);
  root->children.pop_back();

  // set the branch lengths of the children of theRoot
  if (outgrp->isLeaf()) {
    root->length = 0.05 * outgrp->length;
  } else {
    root->length = 0.5 * outgrp->length;
  }
  outgrp->length -= root->length;

  return theRoot;
};

/********************************************************************************
 * @brief option the tree for balance
 *
 * @param root
 * @return Node*
 ********************************************************************************/
void Node::balanceTree(const string &meth, size_t otulvl) {
  // set the depth
  _getDepth();

  // find the root candidates
  vector<Node *> nlist;
  if (otulvl <= taxLevel)
    otulvl = taxLevel + 1;
  findRootCandidates(nlist, otulvl);

  // do the rooting
  if (meth.compare("mdmp") == 0) {
    _mdmpTree(nlist);
  } else if (meth.compare("pmr") == 0) {
    _pmrTree(nlist);
  } else if (meth.compare("mad") == 0) {
    _madTree(nlist);
  } else {
    theInfo("Unkown rooting method: " + meth + ", use MAD instead");
    _madTree(nlist);
  }
};

void Node::findRootCandidates(vector<Node *> &nlist, size_t otulvl) {
  // set the front as the outgroup
  if (children.front()->taxLevel != taxLevel) {
    Node *tmp = children.back();
    children.back() = children.front();
    children.front() = tmp;
  }

  // get the candidate for outgroup
  nlist.emplace_back(children.back());
  children.front()->_findRootCandidates(nlist, otulvl);
  children.back()->_findRootCandidates(nlist, otulvl);
  theInfo("There are " + to_string(nlist.size()) + " root candidates");
}

void Node::_findRootCandidates(vector<Node *> &nlist, size_t otulvl) {
  if (taxLevel < otulvl) {
    for (auto nd : children) {
      nlist.emplace_back(nd);
      if (!nd->isLeaf() && nd->nleaf>1) {
        nd->_findRootCandidates(nlist, otulvl);
      } else {
        nd->otu = true;
      }
    }
  } else {
    otu = true;
  }
}

void Node::_rearrangeOutgroup(Node *np) {
  // check the outgroup if the root
  if (np->parent == NULL)
    return;

  // get the out group
  Node *outgrp = np;
  Node *ogSib = outgrp->_sibling();

  // check the outgroup is already the outgroup
  if (np->parent->parent == NULL) {
    outgrp->length += ogSib->length;
    ogSib->length = 0;
    return;
  }

  // relative nodes: from parent to the origal root
  vector<Node *> nlist;
  do {
    nlist.emplace_back(np->parent);
    np = np->parent;
  } while (np->parent->parent != NULL);

  // get the subroot
  Node *subRoot = nlist.back();
  nlist.pop_back();
  if (subRoot != children.front()) {
    children.front()->swap(children.back());
    subRoot = children.front();
  }

  // set the length of rest and remove it from root
  Node *rest = subRoot->_sibling();
  rest->length += subRoot->length;
  subRoot->length = 0;
  deleteChild(rest);
  // add outgroup to root
  addChild(outgrp);

  // reset the subroot
  if (!nlist.empty()) {
    // clear children and add sibling as a child
    for (auto nd : nlist) {
      nd->children.clear();
      nd->addChild(nd->_sibling());
    }

    // add parent as a child
    for (auto iter = nlist.rbegin() + 1; iter != nlist.rend(); ++iter) {
      (*iter)->addChild((*iter)->parent);
    }
    nlist.back()->addChild(rest);

    // update the status of serial nodes
    for (auto iter = nlist.rbegin(); iter != nlist.rend(); ++iter) {
      (*iter)->_setOneBranch();
      (*iter)->_getDepth();
    }

    // set the top brache as rest
    rest = nlist.front();
  }

  // reset the subroot
  subRoot->children.clear();
  subRoot->addChild(ogSib);
  subRoot->addChild(rest);
  subRoot->_setOneBranch();
  subRoot->_getDepth();
};

Node *Node::_sibling() {
  return parent->children.front() == this ? parent->children.back()
                                          : parent->children.front();
};

/********************************************************************************
 * @brief get minimal ancestor deviation tree
 *
 * @param nlist
 ********************************************************************************/
void Node::_madTree(const vector<Node *> &nlist) {
  // find the minimal depth
  pair<Node *, pair<double, double>> minMAD{
      NULL, make_pair(0, numeric_limits<double>::max())};
  for (auto nd : nlist) {
    // rearrange the tree by set nd as the outgroup
    _rearrangeOutgroup(nd);
    pair<double, double> mad = getMAD();

    // get max pmr
    if (mad.second < minMAD.second.second) {
      minMAD = make_pair(nd, mad);
    }
  }

  // select the best tree
  _rearrangeOutgroup(minMAD.first);
  children.back()->length -= minMAD.second.first;
  children.front()->length = minMAD.second.first;

  theInfo("Rooting Tree by minimal ancestor deviation: " +
          to_string(minMAD.second.second));
}

pair<double, double> Node::getMAD() {
  // get depthes of all otu
  double mad(0);
  vector<double> backdep;
  children.back()->_getOTUad(backdep, mad);
  vector<double> frontdep;
  children.front()->_getOTUad(frontdep, mad);

  // obtian the break point
  double sumx(0);
  double sumy(0);
  for (auto &nda : frontdep) {
    for (auto &ndb : backdep) {
      double dbc = 1 / (nda + ndb);
      double ddbc = dbc * dbc;
      sumx += (ndb - nda) * ddbc;
      sumy += ddbc;
    }
  }
  double doi = sumx * 0.5 / sumy;

  double doj2;
  if (doi < 0) {
    doi = 0;
    doj2 = children.back()->length * 2;
  } else if (doi > children.back()->length) {
    doi = children.back()->length;
    doj2 = 0;
  } else {
    doj2 = (children.back()->length - doi) * 2;
  }

  // add the sanddle ancestor deviation
  for (auto &nda : frontdep) {
    for (auto &ndb : backdep) {
      double dev = (ndb - nda - doj2) / (nda + ndb);
      mad += dev * dev;
    }
  }

  return make_pair(doi, mad);
}

void Node::_getOTUad(vector<double> &theDepth, double &mad) {
  if (otu) {
    theDepth.emplace_back(depth + length);
  } else {
    size_t b = theDepth.size();
    children.front()->_getOTUad(theDepth, mad);
    size_t m = theDepth.size();
    children.back()->_getOTUad(theDepth, mad);
    size_t e = theDepth.size();

    // add the ancestor deviation between intro OTU
    for (size_t i = b; i < m; ++i) {
      for (size_t j = m; j < e; ++j) {
        double rbc = (theDepth[i] - theDepth[j]) / (theDepth[i] + theDepth[j]);
        mad += rbc * rbc;
      }
    }

    // update depth of the descendant OTU
    for (size_t i = b; i < e; ++i) {
      theDepth[i] += length;
    }
  }
}

/********************************************************************************
 * @brief get the maximal relative pairwise midpoint root tree
 *
 * @param nlist
 ********************************************************************************/
void Node::_pmrTree(const vector<Node *> &nlist) {
  // find the minimal depth
  pair<Node *, double> maxPMR{NULL, numeric_limits<double>::min()};
  for (auto nd : nlist) {
    // rearrange the tree by set nd as the outgroup
    _rearrangeOutgroup(nd);
    double pmr = getPMR();

    // get max pmr
    if (pmr > maxPMR.second) {
      maxPMR = make_pair(nd, pmr);
    }
  }

  // select the best tree
  _rearrangeOutgroup(maxPMR.first);
  _setLengthByMidpoint();

  stringstream buf;
  buf << fixed << setprecision(1) << maxPMR.second * 100 << "%";
  theInfo("Rooting Tree by pairwise midpoint root with percent: " + buf.str());
}

double Node::getPMR() {
  // get depthes of all otu
  vector<double> backdep;
  children.back()->_getOTUdepth(backdep);
  vector<double> frontdep;
  children.front()->_getOTUdepth(frontdep);

  // obtian the precent of PMR
  int nPMR(0);
  for (auto &nda : frontdep) {
    for (auto &ndb : backdep) {
      double dio = 0.5 * (ndb - nda);
      if (dio > 0 && dio < children.back()->length) {
        ++nPMR;
      }
    }
  }

  return double(nPMR) / (frontdep.size() * backdep.size());
};

void Node::_getOTUdepth(vector<double> &theDepth) {
  if (otu) {
    theDepth.emplace_back(depth + length);
  } else {
    size_t i = theDepth.size();
    for (auto nd : children) {
      nd->_getOTUdepth(theDepth);
    }

    for (; i < theDepth.size(); ++i) {
      theDepth[i] += length;
    }
  }
}

/********************************************************************************
 * @brief get the minimal depth tree with midpoint and positive length
 *
 * @param nlist
 ********************************************************************************/
void Node::_mdmpTree(const vector<Node *> &nlist) {
  // find the minimal depth
  pair<Node *, double> mindep{NULL, numeric_limits<double>::max()};
  pair<Node *, double> mindepPlus{NULL, numeric_limits<double>::max()};
  for (auto nd : nlist) {
    // rearrange the tree by set nd as the outgroup
    _rearrangeOutgroup(nd);
    // reset the branch
    _setLengthByMidpoint();

    // get mindep
    if (depth < mindep.second) {
      mindep = make_pair(nd, depth);
    }

    // get mindepPlus
    if (depth < mindepPlus.second && children.front()->length > 0 &&
        children.back()->length > 0) {
      mindepPlus = make_pair(nd, depth);
    }
  }

  // select the best tree
  if (mindepPlus.first == NULL) {
    _rearrangeOutgroup(mindep.first);
  } else {
    _rearrangeOutgroup(mindepPlus.first);
  }
  _setLengthByMidpoint();

  theInfo("Rooting Tree by minimal depth: " + to_string(depth));
}

void Node::_setLengthByMidpoint() {
  // points to two branches
  Node *deep = children.front();
  Node *shallow = children.back();
  if (deep->depth < shallow->depth) {
    deep = children.back();
    shallow = children.front();
  }

  // set new branch length (bigger than 0)
  double lengthTotal = deep->length + shallow->length;
  deep->length = 0.5 * (shallow->depth + lengthTotal - deep->depth);
  if (deep->length < 0) {
    deep->length = 0;
    shallow->length = lengthTotal;
  } else {
    shallow->length = lengthTotal - deep->length;
  }

  // update the depth
  _getDepth();
}

void Node::_getDepth() {
  depth = 0;
  if (!isLeaf()) {
    for (auto &nd : children) {
      if (nd->isLeaf()) {
        nd->depth = 0;
        depth += nd->length;
      } else {
        if (std::isnan(nd->depth))
          nd->_getDepth();
        depth += ((nd->length + nd->depth) * (nd->nleaf + nd->nxleaf));
      }
    }
    depth /= double(nleaf+nxleaf);
  }
};

/********************************************************************************
 * @brief reset rearrange the tree by change the outgroup for unroot tree
 *
 * @param str
 * @return Node*
 ********************************************************************************/

Node *Node::resetroot(const string &str) {
  vector<Node *> leafs;
  getLeafs(leafs);
  Node *outgrp = NULL;
  for (Node *nd : leafs) {
    if (lastNameNoRank(nd->name) == str) {
      outgrp = nd;
      break;
    }
  }

  if (outgrp == NULL)
    return outgrp;
  return resetroot(outgrp);
}

Node *Node::resetroot(Node *np) {
  Node *outgrp = np;
  vector<Node *> nlist;
  // relative nodes: from parent to the origal root
  while ((*np).parent != NULL) {
    nlist.emplace_back((*np).parent);
    np = (*np).parent;
  }

  // reset the length and toplogy
  for (auto iter = nlist.rbegin(); iter != nlist.rend() - 1; ++iter) {
    auto next = iter + 1;
    (**iter).length = (**next).length;
    (**next).length = NAN;
    (**iter).deleteChild(*next);
    (**next).parent = NULL;
    (**next).addChild(*iter);
  }
  Node *theRoot = nlist.front();

  // set the outgroup as the last child
  auto iter = (*theRoot).children.rbegin();
  if ((*iter) != outgrp) {
    for (++iter; iter != (*theRoot).children.rend(); ++iter) {
      if ((*iter) == outgrp) {
        *iter = (*theRoot).children.back();
        (*theRoot).children.back() = outgrp;
      }
    }
  }

  return theRoot;
};

/********************************************************************************
 * @brief section for option on nwk file
 *
 * @param file
 ********************************************************************************/
function<string(Node *)> Node::nwkname = [](Node *nd) {
  return nd->name.substr(nd->name.find_first_of('|') + 1);
};

void Node::_outnwk(ostream &os) {
  if (!isLeaf()) {
    os << "(";
    auto iter = children.begin();
    (*iter)->_outnwk(os);
    for (++iter; iter != children.end(); ++iter) {
      os << ",";
      (*iter)->_outnwk(os);
    }
    os << ")";
  }

  // output the name of node
  string namestr = nwkname(this);
  os << ((namestr.find(' ') == string::npos) ? namestr : ('"' + namestr + '"'));

  // output the length of node
  if (!std::isnan(length))
    os << ":" << fixed << setprecision(5) << length;
};

void Node::outnwk(const string &file) {
  ofstream os(file);
  if (!os) {
    cerr << "Open " << file << " for write failed" << endl;
    exit(3);
  }

  outnwk(os);
  os.close();
};

void Node::outitol(const string &file) {
  ofstream os(file);
  if (!os) {
    cerr << "Open " << file << " for write failed" << endl;
    exit(3);
  }

  // output the taxid
  int theId = 0;
  updateId(theId);
  nwkname = [](Node *nd) {
    if (nd->isLeaf()) {
      return to_string(nd->id);
    } else {
      return "I" + to_string(nd->id);
    }
  };

  outnwk(os);
  os.close();
};

void Node::outnwk(ostream &os) {
  _outnwk(os);
  os << ";" << endl;
};

void Node::_nwkItem(const string &str) {
  vector<string> words;
  separateWord(words, str, ":");
  if (isLeaf()) {
    name = words.front();
  } else {
    try {
      bootstrap = str2double(words.front());
    } catch (exception e) {
      name = words.front();
    }
  }
  length = str2double(words.back());
};

void Node::_innwk(istream &is) {
  Node *np = new Node;
  addChild(np);
  string brStr;

  while (is.good()) {
    char c = is.get();
    if (c == ',') {
      if (!brStr.empty()) {
        np->_nwkItem(brStr);
        brStr.clear();
      }

      np = new Node;
      addChild(np);
    } else if (c == ')') {
      if (!brStr.empty()) {
        np->_nwkItem(brStr);
        brStr.clear();
      }
      break;
    } else if (c == '(') {
      np->_innwk(is);
    } else if (c != '"' && c != '\n' && c != '\t' && c != '\'') {
      brStr.push_back(c);
    }
  }
};

void Node::innwk(istream &is) {

  char c = is.get();
  string brStr;
  while (c != ';') {
    if (c == '(') {
      _innwk(is);
    } else if (c != '"' && c != '\n' && c != '\t' && c != '\'') {
      brStr.push_back(c);
    }

    if (is.good()) {
      c = is.get();
    } else {
      cerr << "some wrong in the nwk file" << endl;
      exit(1);
    }
  }
};

void Node::innwk(const string &file) {
  ifstream inwk(file);
  if (!inwk) {
    cerr << "Cannot found the input newick file " << file << endl;
    exit(4);
  }

  innwk(inwk);
  inwk.close();
}

void Node::checkUnclassified() {
  if (isLeaf()) {
    if (unclassified) {
      nxleaf = 1;
    } else {
      nleaf = 1;
    }
  } else {
    nxleaf = 0;
    nleaf = 0;
    for (const auto &nd : children) {
      nd->checkUnclassified();
      nxleaf += nd->nxleaf;
      nleaf += nd->nleaf;
    }
  }

  if (nleaf == 0) {
    unclassified = true;
  }
};

void Node::checkUploaded() {
  if (isLeaf()) {
    if (name.find(".UPLOAD") != string::npos) {
      nUpLeaf = 1;
      uploaded = true;
    }
  } else {
    bool isUpload = true;
    for (const auto &nd : children) {
      (*nd).checkUploaded();
      nUpLeaf += (*nd).nUpLeaf;
      isUpload = isUpload && (*nd).uploaded;
    }
    uploaded = isUpload;
  }
};

void Node::getDefineLeafs(vector<Node *> &nodes) {
  if (isLeaf()) {
    if (!unclassified)
      nodes.emplace_back(this);
  } else {
    for (const auto &nd : children)
      (*nd).getDefineLeafs(nodes);
  }
};

/********************************************************************************
 * @brief option on json file
 *
 * @param file
 ********************************************************************************/

void Node::injson(const string &file) {
  ifstream ijson(file.c_str());
  if (!ijson) {
    cerr << " Cannot found the input file " << file << endl;
    exit(4);
  }
  injson(ijson);
}

void Node::injson(istream &is) {
  is.get();
  _injson(is);
};

void Node::_injson(istream &is) {
  char c = is.get();
  bool isValue(false);
  string key, value;
  while (c != '}') {
    if (c == '"') {
      if (isValue)
        _getStr(is, value);
      else
        _getStr(is, key);
    } else if (c == ':') {
      isValue = true;
    } else if (c == '{') {
      Node *nd = new Node;
      addChild(nd);
      (*nd)._injson(is);
    } else if (c == ',' && isValue) {
      _getKeyValue(key, value);
      isValue = false;
    }
    c = is.get();
  }
  if (isValue)
    _getKeyValue(key, value);
};

void Node::_getKeyValue(string &key, string &value) {
  if (key.compare("name") == 0) {

    // find the arch
    size_t lpos = value.find('{');
    size_t rpos = value.find('}');
    size_t spos = value.find('/');
    size_t ppos = value.find('+');

    // get the name
    name = value.substr(0, lpos);

    // get the unclassified number
    if (ppos == string::npos) {
      nxleaf = 0;
    } else {
      nxleaf = str2int(value.substr(ppos + 1, rpos - ppos - 1));
      rpos = ppos;
    }

    // get the classified number and taxonomy size
    if (spos == string::npos) {
      taxSize = str2int(value.substr(lpos + 1, rpos - lpos - 1));
      nleaf = taxSize;
    } else {
      nleaf = str2int(value.substr(lpos + 1, spos - lpos - 1));
      taxSize = str2int(value.substr(spos + 1, rpos - spos - 1));
    }
  } else if (key.compare("length") == 0) {
    length = str2double(value);
  } else if (key.compare("upload") == 0) {
    uploaded = true;
  } else if (key.compare("unclassified")) {
    unclassified = true;
  }
};

void Node::_getStr(istream &is, string &str) {
  str.clear();
  char c = is.get();
  while (c != '"') {
    str += c;
    c = is.get();
  }
};

void Node::outjson(const string &file) {
  ofstream os(file);
  if (!os) {
    cerr << "Open " << file << " for write failed" << endl;
    exit(3);
  }

  outjson(os);
  os.close();
};

void Node::outjson(ostream &os) {

  os << "{";
  if (!name.empty()) {
    os << '"' << "name"
       << "\":\"" << name << "{" << nleaf;
    if (nleaf > 0 && nleaf < taxSize)
      os << "/" << taxSize;
    if (nxleaf != 0) {
      os << "+" << nxleaf;
      // if(nleaf==0 && nxleaf < taxSize)
      // 	os << "/" << taxSize;
    }
    os << "}\"";
  }

  if (uploaded)
    os << ',' << '"' << "upload"
       << "\":\""
       << "true" << '"';

  if (unclassified)
    os << ',' << '"' << "unclassified"
       << "\":\""
       << "true" << '"';

  if (nleaf == taxSize) {
    os << ',' << '"' << "ntype"
       << "\":\""
       << "Coincide" << '"';
  } else if (nleaf == 0 && nxleaf == taxSize) {
    os << ',' << '"' << "ntype"
       << "\":\""
       << "Coincide" << '"';
  }

  if (!std::isnan(length))
    os << ',' << '"' << "length"
       << "\":\"" << fixed << setprecision(5) << length << '"';

  if (!isLeaf()) {
    os << ",\"children\":[";

    Children::const_iterator iter = children.begin();
    (*(*iter)).outjson(os);
    for (++iter; iter != children.end(); ++iter) {
      os << ",";
      (*(*iter)).outjson(os);
    }
    os << "]";
  }
  os << "}";
};

void Node::outjsonAbbr(ostream &os) {

  os << "{";
  if (!name.empty()) {
    os << '"' << "n"
       << "\":\"" << name << "{" << nleaf;
    if (nleaf < taxSize)
      os << "/" << taxSize;
    if (nxleaf != 0) {
      os << "+" << nxleaf;
      // if(nleaf==0 && nxleaf < taxSize)
      // 	os << "/" << taxSize;
    }
    os << "}\"";
  }

  if (uploaded)
    os << ',' << '"' << "u"
       << "\":\""
       << "1" << '"';

  if (unclassified)
    os << ',' << '"' << "x"
       << "\":\""
       << "1" << '"';

  if (nleaf == taxSize)
    os << ',' << '"' << "t"
       << "\":\""
       << "C" << '"';
  else if (nleaf == 0 && nxleaf == taxSize) {
    os << ',' << '"' << "t"
       << "\":\""
       << "C" << '"';
  }

  // if(!isnan(length))
  // 	os << ',' << '"' << "l" << "\":\"" << fixed << setprecision(5) << length
  // <<'"';

  if (!isLeaf()) {
    os << ",\"c\":[";

    Children::const_iterator iter = children.begin();
    (*(*iter)).outjsonAbbr(os);
    for (++iter; iter != children.end(); ++iter) {
      os << ",";
      (*(*iter)).outjsonAbbr(os);
    }
    os << "]";
  }
  os << "}";
};

void Node::renewId(const unordered_map<string, size_t> &mgi) {
  if (isLeaf()) {
    id = mgi.find(name)->second;
  } else {
    for (auto &nd : children)
      (*nd).renewId(mgi);
  }
};

void Node::updateId(int &theId) {
  id = ++theId;
  if (!isLeaf()) {
    for (auto &nd : children)
      (*nd).updateId(theId);
  }
};

/********************************************************************************
 * @brief
 *
 * @param file
 ********************************************************************************/
void Node::outPrediction(const string &file) {
  ofstream os(file);
  if (!os) {
    cerr << "Open " << file << " for write failed" << endl;
    exit(3);
  }

  outPrediction(os);
  os.close();
};

void Node::outPrediction(ostream &os) {
  vector<Node *> leafs;
  getLeafs(leafs);
  for (auto &nd : leafs) {
    if ((*nd).unclassified) {
      string p;
      (*nd)._getPrediction(p);
      p.erase(remove(p.begin(), p.end(), '|'), p.end());
      os << lastNameNoRank(nd->name) << "\t" << p << endl;
    }
  }
};

void Node::_getPrediction(string &p) {
  if (parent == NULL) {
    p = name;
  } else {
    if ((*parent).nleaf == (*parent).taxSize) {
      p = (*parent).name;
    } else {
      (*parent)._getPrediction(p);
    }
  }
};

void Node::reinitTree() {
  if (isLeaf()) {
    name.erase(name.find('|'), 1);
  } else {
    name = "";
    for (const auto &nd : children)
      (*nd).reinitTree();
  }

  nxleaf = 0;
  nleaf = 0;
  taxSize = 0;
  unclassified = false;
}

void Node::getUndefineLeafs(vector<Node *> &nodes) {
  if (children.empty()) {
    if (unclassified)
      nodes.emplace_back(this);
  } else {
    for (auto nd : children)
      (*nd).getUndefineLeafs(nodes);
  }
};

void Node::getUndefineNames(vector<string> &names) {
  if (children.empty()) {
    if (unclassified)
      names.emplace_back(name);
  } else {
    for (auto nd : children)
      (*nd).getUndefineNames(names);
  }
};

/********************************************************************************
 * @brief set the lineage and annotate taxonomy of tree
 *
 ********************************************************************************/
void Node::setOneLeaf(const string &nm, bool def) {
  if (def) {
    nleaf = 1;
    unclassified = false;
  } else {
    nxleaf = 1;
    unclassified = true;
  }
  name = nm;
  taxLevel = nRanks(nm);
};

void Node::_setOneBranch() {

  nxleaf = 0;
  nleaf = 0;
  vector<Node *> defNode;
  vector<Node *> undefNode;
  for (auto &nd : children) {
    nxleaf += nd->nxleaf;
    nleaf += nd->nleaf;
    if (nd->unclassified) {
      undefNode.emplace_back(nd);
    } else {
      defNode.emplace_back(nd);
    }
  }

  // set unclassified and name
  if (nleaf == 0) {
    // for unclassified node
    unclassified = true;
    name = _getBranchName(undefNode);
  } else {
    // for defined node
    name = _getBranchName(defNode);
  }

  // set tax level
  taxLevel = nRanks(name);
}

void Node::setAllBranches() {
  if (!isLeaf()) {
    for (auto &child : children)
      child->setAllBranches();
    _setOneBranch();
  }
}

void Node::setBranchLineage() {
  if (!isLeaf()) {
    vector<Node *> defNode;
    vector<Node *> undefNode;
    for (auto &child : children) {
      child->setBranchLineage();

      if (child->unclassified) {
        undefNode.emplace_back(child);
      } else {
        defNode.emplace_back(child);
      }
    }

    if (defNode.size() > 0) {
      name = _getBranchName(defNode);
    } else {
      name = _getBranchName(undefNode);
    }
  }

  // set tax level
  taxLevel = nRanks(name);
};

string Node::_getBranchName(const vector<Node *> &nds) {
  string nodeName;
  if (nds.size() == 1) {
    nodeName = delStrain(nds.front()->name);
  } else {
    auto iter = nds.begin();
    nodeName = (*iter++)->name;
    do {
      string lngB = (*iter++)->name;
      nodeName = commonLineage(nodeName, lngB);
    } while (iter != nds.end());
  }

  return nodeName;
};

int Node::nClade() {
  if (parent == NULL) {
    return taxLevel;
  } else {
    size_t pTaxLevel = parent->taxLevel;
    if (unclassified && !parent->unclassified)
      pTaxLevel = 0;
    return taxLevel - pTaxLevel;
  }
};
