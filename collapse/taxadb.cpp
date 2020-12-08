/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-12-05 15:07:07
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-07 20:54:00
 */

#include "taxadb.h"
/********************************************************************************
 * @brief for a node of taxonomy system
 *
 ********************************************************************************/
ostream &operator<<(ostream &os, const TaxonNode &tax) {
  os << tax.taxid << "\t" << tax.level << "\t" << tax.name << "\t"
     << tax.lineage;
  return os;
};

/********************************************************************************
 * @brief Construct a new TaxaDB:: TaxaDB object
 *
 * @param fnode path of nodes.dmp
 ********************************************************************************/
TaxaDB::TaxaDB(const string &path) {
  string dbpath = path;
  addsuffix(dbpath, '/');
  string fnode = dbpath + "nodes.dmp";
  string fname = dbpath + "names.dmp";

  _readNodeDump(fnode);
  _readNameDump(fname);

  theInfo("INITIAL Section: initial database by files in " + dbpath);
};

void TaxaDB::_readNodeDump(const string &file) {
  ifstream ifs(file);
  if (!ifs) {
    cerr << "cannot open file " << file << endl;
    exit(1);
  }

  for (string line; getline(ifs, line);) {
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
  ifs.close();
};

void TaxaDB::_readNameDump(const string &file) {
  ifstream ifs(file);
  if (!ifs) {
    cerr << "cannot open file " << file << endl;
    exit(1);
  }

  for (string line; getline(ifs, line);) {
    line = trim(line);
    if (line.empty())
      continue;

    vector<string> words;
    separateWord(words, line, "|");
    size_t tid = stoi(words[0]);
    string nm = trim(words[1]);
    string type = trim(words[3]);
    if (type == "scientific name") {
      nm = regex_replace(nm, nonWord, "_");
      if (nm.back() == '_')
        nm.pop_back();
      taxlist.at(tid).name = nm;
      taxlist.at(tid).lowerName = toLower(nm);
    }
  }

  ifs.close();
};

// set rank name map
void TaxaDB::resetRankMap(const string &file) { hdl.setRankByFile(file); };

/********************************************************************************
 * @brief search the lineage of the query name
 *
 * @param qry name of query
 * @return string return the lineage string
 ********************************************************************************/
string TaxaDB::search(string qry) {

  // format the query
  qry = regex_replace(qry, nonWord, "_");
  if (qry.back() == '_')
    qry.pop_back();
  qry = toLower(qry);

  // query the lineage by Taxid
  smatch matchs;
  regex_search(qry, matchs, regTaxID);
  if (!matchs.empty()) {
    TaxonNode tax = getLineage(stoi(matchs[1].str()));
    if (tax.taxid != 0) {
      return tax.lineage;
    }
  }

  // query the lineage by subname  
  size_t pos;
  do {
    TaxonNode tax = getLineage(qry);
    if (tax.taxid != 0) {
      return tax.lineage;
    }
    pos = qry.rfind("_");
    qry = qry.substr(0, pos);
  } while (pos != string::npos);

  theInfo("Cannot find the lineage for " + qry);
  return "";
};

// get lineage by taxid
TaxonNode TaxaDB::getLineage(size_t id) {
  _getlng(id);
  return taxlist[id];
};

// get lineage by name
TaxonNode TaxaDB::getLineage(const string &nm) {
  for (auto &tax : taxlist) {
    if (tax.lowerName == nm) {
      _getlng(tax.taxid);
      return tax;
    }
  }
  return taxlist.front();
};

// export all lineages to a ostream
void TaxaDB::exportLineage(ostream &os) {
  for (auto &tax : taxlist) {
    if (tax.taxid != 0) {
      _getlng(tax.taxid);
      os << tax << "\n";
    }
  }
};

// recursively get the lineage string
string TaxaDB::_getlng(size_t tid) {
  auto &tax = taxlist[tid];
  if (tax.taxid == tax.parent) {
    return "";
  } else if (tax.lineage.empty()) {
    tax.lineage =
        _getlng(tax.parent) + hdl.lineageString(make_pair(tax.level, tax.name));
  }
  return tax.lineage;
};
