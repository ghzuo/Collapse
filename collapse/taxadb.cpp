/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-12-05 15:07:07
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-06 17:11:08
 */

#include "taxadb.h"

ostream &operator<<(ostream &os, const TaxonNode &tax) {
  os << tax.taxid << "\t" << tax.level << "\t" << tax.name << "\t"
     << tax.lineage;
  return os;
};

TaxaDB::TaxaDB(const string &fnode, const string &fname) {
  _readNodeDump(fnode);
  _readNameDump(fname);

  theInfo("INITIAL Section: initial database by " + fnode + " and " + fname);
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
    }
  }

  ifs.close();
};

// set rank name map
void TaxaDB::resetRankMap(const string &file) { hdl.setRankByFile(file); };

// for query
string TaxaDB::search(const string &qry) {

  // directly search the name
  TaxonNode tax = getLineage(qry);
  if (tax.taxid != 0) {
    return tax.lineage;
  }

  // query the lineage by Taxid
  smatch matchs;
  regex_search(qry, matchs, regTaxID);
  if (!matchs.empty()) {
    TaxonNode tax = getLineage(stoi(matchs[1].str()));
    if (tax.taxid != 0) {
      return tax.lineage;
    }
  }

  // query the lineage by specie name
  regex_search(qry, matchs, regSpeciesName);
  if (!matchs.empty()) {
    TaxonNode tax = getLineage(matchs[1].str());
    if (tax.taxid != 0) {
      return tax.lineage;
    }
  }

  // query the lineage by other taxon name
  regex_search(qry, matchs, regGenusName);
  if (!matchs.empty()) {
    TaxonNode tax = getLineage(matchs[1].str());
    if (tax.taxid != 0) {
      return tax.lineage;
    }
  }

  theInfo("Cannot find the lineage for " + qry);
  return "";
};

// get lineage
TaxonNode TaxaDB::getLineage(size_t id) {
  _getlng(id);
  return taxlist[id];
};

TaxonNode TaxaDB::getLineage(const string &nm) {
  for (auto &tax : taxlist) {
    if (tax.name == nm) {
      _getlng(tax.taxid);
      return tax;
    }
  }
  theInfo("Cannot find a taxon named as " + nm);
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
