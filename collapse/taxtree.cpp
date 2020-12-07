/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-03-17 15:39:23
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-05 21:58:16
 */

#include "taxtree.h"
const size_t N_FORKS(2);

/*******************************************************************/
/********* Member Functions For Node Class *************************/
/*******************************************************************/
Node::Node()
    : name(""), id(0), length(NAN), parent(NULL), taxSize(0), taxLevel(0),
      nleaf(0), nxleaf(0), unclassified(false), uploaded(false) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n)
    : name(""), id(n), length(NAN), parent(NULL), taxSize(0), taxLevel(0),
      nleaf(0), nxleaf(0), unclassified(false), uploaded(false) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n, const string &str)
    : name(str), id(n), length(NAN), parent(NULL), taxSize(0), taxLevel(0),
      nleaf(0), nxleaf(0), unclassified(false), uploaded(false) {
  children.reserve(N_FORKS);
};

Node::Node(size_t n, const vector<Node *> &vn)
    : name(""), id(n), length(NAN), parent(NULL), children(vn), taxSize(0),
      nleaf(0), nxleaf(0), unclassified(false), uploaded(false) {
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

/********************************************************************************
 * @brief tools for set root of the unroot tree
 *
 * @param str
 * @return Node*
 ********************************************************************************/
void Node::updateRootedTree() {
  theInfo("This tree is a rooted treee, keep as it is");

  // set branch lineage: type, fullname,  taxLevel
  checkUnclassified();
  setBranchLineage();
}

Node *Node::rootingDirect(){
  Node* theTree = _forceRooting(this);
  theTree->checkUnclassified();
  theTree->setBranchLineage();
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
  theTree->checkUnclassified();
  theTree->setBranchLineage();
  return theTree;
};

Node *Node::rootingByTaxa() {

  // preset the branch attraction
  Node *theTree = this;
  theTree->checkUnclassified();
  theTree->setBranchLineage();

  // find the candidate of outgroup
  vector<Node *> clades;
  for (auto &nd : theTree->children) {
    nd->_getOutgrpCandidates(clades);
  }

  // root the tree by the highest rank clade and length
  if (!clades.empty()) {
    Node *outgrp(clades.front());

    for (auto &nd : clades) {
      if (outgrp->taxLevel > nd->taxLevel ||
          (outgrp->taxLevel == nd->taxLevel && outgrp->length < nd->length)) {
        outgrp = nd;
      }
    }

    // record the change branch
    Node *chgBranch = outgrp->parent->parent;
    // rearrange the tree (still unroot tree)
    theTree = resetroot(outgrp);
    // renew the update branch
    chgBranch->checkUnclassified();
    chgBranch->setBranchLineage();

    // find a good outgroup (the last item of children)
    // by higest rank of common lineage or longest branch length
    size_t minComLev = nRanks(commonLineage(chgBranch->name, outgrp->name));
    size_t maxLength = outgrp->length;
    auto outIter = children.rbegin();
    auto chgIter = children.rend();
    for (auto iter = theTree->children.rbegin() + 1;
         iter != theTree->children.rend(); ++iter) {
      if (*iter != chgBranch) {
        size_t comLev = nRanks(commonLineage(chgBranch->name, (*iter)->name));
        if ((comLev < minComLev) ||
            ((comLev == minComLev) && ((*iter)->length > maxLength))) {
          outIter = iter;
          minComLev = comLev;
          maxLength = (*iter)->length;
        }
      } else {
        chgIter = iter;
      }
    }

    // swap new outgroup to the back
    if ((*chgIter)->length > maxLength) {
      auto tmp = *chgIter;
      *chgIter = outgrp;
      theTree->children.back() = tmp;
    } else if (outgrp != *outIter) {
      auto tmp = *outIter;
      *outIter = outgrp;
      theTree->children.back() = tmp;
    }
  }

  // rooting the unroot tree to a rooted tree
  theTree = _forceRooting(theTree);

  // renew the added node and the root
  theTree->children.back()->_renewNode();
  theTree->_renewNode();

  return theTree;
};

void Node::_renewNode() {
  nxleaf = 0;
  nleaf = 0;
  for (auto &nd : children) {
    nxleaf += nd->nxleaf;
    nleaf += nd->nleaf;
  }
  if (nleaf == 0)
    unclassified = true;
}

bool Node::_getOutgrpCandidates(vector<Node *> &clades) {
  if (unclassified) {
    return false;
  } else if (isLeaf() || nClade() > 0) {
    return true;
  } else {
    // the candidates are the nodes whose all siblings are clades
    bool isCandidate(true);
    for (auto nd : children) {
      isCandidate = nd->_getOutgrpCandidates(clades) && isCandidate;
    }
    if (isCandidate) {
      clades.insert(clades.end(), children.begin(), children.end());
    }
    return false;
  }
};

Node *Node::_forceRooting(Node *root) {

  // testing whether a rooted tree
  if(root->children.size()==2){
    theInfo("This is a rooted tree");
    return root;
  }

  // add a super root
  Node *theRoot = new Node();

  // add the last child of node as the outgroup of theRoot
  Node *outgrp = root->children.back();
  // the branch length are divided equally to two branch
  outgrp->length *= 0.5;
  theRoot->children.push_back(outgrp); // outgrp at front
  root->children.pop_back();

  // the origal node as another child of the Root
  theRoot->children.push_back(root);
  root->length = outgrp->length;

  return theRoot;
};

Node *Node::resetroot(const string &str) {
  vector<Node *> leafs;
  getLeafs(leafs);
  Node *outgrp = NULL;
  for (Node *nd : leafs) {
    if (_lastName(nd->name) == str) {
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
  string forename = name.substr(name.find_first_of('|') + 1);
  os << ((name.find(' ') == string::npos) ? forename : ('"' + forename + '"'));

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

/********************************************************************************
 * @brief
 *
 * @param file
 ********************************************************************************/
string Node::getStrainName() {
  return name.substr(name.find_last_of('>') + 1);
};

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
      os << (*nd).getStrainName() << "\t" << p << endl;
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
  taxLevel = 0;
  for (auto c : name) {
    if (c == '<')
      taxLevel++;
  }
};

string Node::_getBranchName(const vector<Node *> &nds) {
  string nodeName;
  if (nds.size() == 1) {
    nodeName = nds.front()->name;
    size_t npos = nodeName.find("<T>");
    if (npos != string::npos) {
      nodeName = nodeName.substr(0, npos);
    }
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
    return 1;
  } else {
    size_t pTaxLevel = parent->taxLevel;
    if (unclassified && !parent->unclassified)
      pTaxLevel = 0;
    return taxLevel - pTaxLevel;
  }
};
