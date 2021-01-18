/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-12-05 15:07:07
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-12 13:47:24
 */

#include "taxadb.h"
/********************************************************************************
 * @brief for a node of taxonomy system
 *
 ********************************************************************************/
ostream &operator<<(ostream &os, const TaxonNode &tax) {
  os << tax.taxid << "\t" << tax.level << "\t" << tax.name;
  return os;
};

/********************************************************************************
 * @brief Construct a new TaxaDB:: TaxaDB object
 *
 * @param fnode path of nodes.dmp
 ********************************************************************************/
TaxaDB::TaxaDB(const string &path) {
  // create the TaxaRank object
  rank = TaxaRank::create();

  // set the system
  if (isDirectory(path)) {

    ifstream fnd(addsuffix(path, '/') + "nodes.dmp");
    _readNodeDump(fnd);
    fnd.close();

    ifstream fnm(addsuffix(path, '/') + "names.dmp");
    _readNameDump(fnm);
    fnm.close();

    theInfo("INITIAL Section: initial database by files in directory " + path);
  } else if (hasSuffix(path, ".tar.gz")) {
    tgz4taxdb(path);
    theInfo("INITIAL Section: initial database by file " + path);
  } else {
    readTable(path);
    theInfo("INITIAL Section: initial database by cache file " + path);
  }
};

/********************************************************************************
 * @brief intial function for database by dump files
 *
 ********************************************************************************/
void TaxaDB::tgz4taxdb(const string &fname) {
  gzFile fp;
  if ((fp = gzopen(fname.c_str(), "rb")) == NULL) {
    cerr << "Cannot found: \"" << fname << '"' << endl;
    exit(1);
  }

  map<string, string> files{{"names.dmp", ""}, {"nodes.dmp", ""}};
  size_t nReadfile(0);
  while (nReadfile < files.size()) {
    TarRecord rc;
    gzread(fp, rc.buff, sizeof(rc));
    if (rc.header.name[0] == 0)
      break;
    size_t nsize = oct2size(rc.header.size);
    auto iter = files.find(rc.header.name);

    // read file
    if (iter != files.end()) {
      tgzReadFile(fp, nsize, iter->second);
      nReadfile++;
    } else {
      size_t bsize = ((nsize + RECORDSIZE - 1) / RECORDSIZE) * RECORDSIZE;
      gzseek(fp, bsize, SEEK_CUR);
    }
  }

  istringstream snd(files["nodes.dmp"]);
  _readNodeDump(snd);

  istringstream snm(files["names.dmp"]);
  _readNameDump(snm);
};

void TaxaDB::_readNodeDump(istream &is) {
  for (string line; getline(is, line);) {
    line = trim(line);
    if (line.empty())
      continue;

    vector<string> words;
    separateWord(words, line, "|");
    TaxonNode nd;
    nd.taxid = stoi(words[0]);
    nd.parent = stoi(words[1]);
    nd.level = regex_replace(trim(words[2]), nonWord, "_");
    if (taxlist.size() <= nd.taxid)
      taxlist.resize(nd.taxid + 1);
    taxlist.at(nd.taxid) = nd;
  }
};

void TaxaDB::_readNameDump(istream &is) {
  for (string line; getline(is, line);) {
    line = trim(line);
    if (line.empty())
      continue;

    vector<string> words;
    separateWord(words, line, "|");
    size_t tid = stoi(words[0]);
    string nm = trim(words[1]);
    string type = trim(words[3]);
    if (type == "scientific name") {
      taxlist.at(tid).name = goodname(nm);
    } else {
      taxlist.at(tid).nicks.push_back(idname(nm));
    }
  }
};

string TaxaDB::goodname(const string &str) {
  string nm(str);
  nm = regex_replace(nm, nonWord, "_");
  if (nm.back() == '_')
    nm.pop_back();
  return nm;
}

string TaxaDB::idname(const string &str) {
  string nm(str);
  nm = regex_replace(nm, nonWord, "");
  return toLower(nm);
}

/********************************************************************************
 * @brief read and write the gnuzip table file of the database
 *
 ********************************************************************************/
void TaxaDB::readTable(const string &fname) {
  gzFile fp;
  if ((fp = gzopen(fname.c_str(), "rb")) == NULL) {
    cerr << "Cannot found: \"" << fname << '"' << endl;
    exit(1);
  }

  string line;
  gzline(fp, line);
  taxlist.resize(stoi(line));

  while (gzline(fp, line) != -1) {
    vector<string> words;
    size_t ncol = separateWord(words, line, ";");
    size_t tid = stoi(words[0]);
    auto &tnd = taxlist[tid];
    tnd.taxid = tid;
    tnd.parent = stoi(words[1]);
    tnd.level = words[2];
    tnd.name = words[3];
    if (ncol == 5)
      separateWord(tnd.nicks, words[4], ",");
  }
};

void TaxaDB::writeTable(const string &fname) {

  gzFile fp;
  if ((fp = gzopen(addsuffix(fname, ".gz").c_str(), "wb")) == NULL) {
    cerr << "Error happen on write database: " << fname << endl;
    exit(1);
  }

  // output the data into a string stream at first
  string head = to_string(taxlist.size()) + "\n";
  gzputs(fp, head.c_str());
  for (auto &tax : taxlist) {
    if (tax.taxid != 0) {
      string str;
      str += to_string(tax.taxid);
      str += ";" + to_string(tax.parent);
      str += ";" + tax.level;
      str += ";" + tax.name;
      if (!tax.nicks.empty())
        str += ";" + strjoin(tax.nicks.begin(), tax.nicks.end(), ',');
      str += "\n";
      gzputs(fp, str.c_str());
    }
  }
  gzclose(fp);
};

/********************************************************************************
 * @brief search the lineage of the query name
 *
 * @param qry name of query
 * @return string return the lineage string
 ********************************************************************************/
string TaxaDB::search(size_t tid) {
  TaxonNode tnd = getLineage(tid);
  if (tnd.taxid) {
    theInfo("Cannot find taxid " + to_string(tid) + " in database");
    return "";
  } else {
    return tnd.lineage;
  }
};

string TaxaDB::search(const string &nm) {

  // format the query
  string qry = toLower(goodname(nm));

  // query the lineage by Taxid
  smatch matchs;
  if (regex_search(qry, matchs, regTaxID)) {
    TaxonNode tax = getLineage(stoi(matchs[1].str()));
    if (tax.taxid != 0) {
      return tax.lineage;
    }
  }

  // query the lineage by subname
  vector<string> qlist;
  _candidateQuery(qry, qlist);
  for (auto &q : qlist) {
    TaxonNode tax = getLineage(_name2id(q));
    if (tax.taxid != 0) {
      return tax.lineage;
    }
  }

  // cannot find the lineage of the genome
  theInfo("Cannot find the lineage for " + nm);
  return "";
};

// get lineage by taxid
TaxonNode TaxaDB::getLineage(size_t id) {
  if (id > taxlist.size())
    return taxlist.front();
  _getlng(id);
  return taxlist[id];
};

// get lineage by name
size_t TaxaDB::_name2id(const string &nm) {

  if (nmlist.empty())
    _initNameMap();
  auto iter = nmlist.find(nm);
  if (iter == nmlist.end()) {
    return 0;
  } else {
    return iter->second;
  }
};

size_t TaxaDB::_candidateQuery(const string &nm, vector<string> &qlist) {
  // for the case from the begin
  string qry(nm);
  size_t pos;
  do {
    qlist.push_back(idname(qry));
    pos = qry.rfind("_");
    qry = qry.substr(0, pos);
  } while (pos != string::npos);

  // for the case without first section name
  qry = nm.substr(nm.find("_") + 1);
  do {
    qlist.push_back(idname(qry));
    pos = qry.rfind("_");
    qry = qry.substr(0, pos);
  } while (pos != string::npos);

  return qlist.size();
};

// initial name search map
void TaxaDB::_initNameMap() {
  for (auto &tax : taxlist) {
    nmlist[idname(tax.name)] = tax.taxid;
    for (auto &nm : tax.nicks) {
      nmlist[nm] = tax.taxid;
    }
  }
}

// recursively get the lineage string
string TaxaDB::_getlng(size_t tid) {
  auto &tax = taxlist[tid];
  if (tax.taxid == tax.parent) {
    return "";
  } else if (tax.lineage.empty()) {
    tax.lineage =
        _getlng(tax.parent) + rank->rankString(make_pair(tax.level, tax.name));
  }
  return tax.lineage;
};

/********************************************************************************
 * @brief for search all taxon nodes of lineage
 *
 ********************************************************************************/
void TaxaDB::query(size_t tid, vector<TaxonNode> &tnlist) {
  if (getFullLineage(tid, tnlist)) {
    theInfo("Cannot find the taxid " + to_string(tid) + " in the database");
  }
};

void TaxaDB::query(const string &nm, vector<TaxonNode> &tnlist) {
  // format the query
  string qry = toLower(goodname(nm));

  // query the lineage by Taxid
  smatch matchs;
  if (regex_search(qry, matchs, regTaxID)) {
    if (getFullLineage(stoi(matchs[1].str()), tnlist)) {
      return;
    }
  }

  // query the lineage by subname
  vector<string> qlist;
  _candidateQuery(qry, qlist);
  for (auto &q : qlist) {
    if (getFullLineage(_name2id(q), tnlist)) {
      return;
    }
  }

  // cannot find the lineage of the genome
  theInfo("Cannot find the lineage for " + nm);
};

bool TaxaDB::getFullLineage(size_t tid, vector<TaxonNode> &tnlist) {
  if (tid > taxlist.size())
    return false;
  _reQuery(tid, tnlist);
  if (tnlist.empty())
    return false;
  return true;
};

void TaxaDB::_reQuery(size_t tid, vector<TaxonNode> &tnlist) {
  auto &tax = taxlist[tid];
  if (tax.taxid != tax.parent) {
    _reQuery(tax.parent, tnlist);
    tnlist.emplace_back(tax);
  }
};
