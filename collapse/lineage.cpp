/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-09-01 13:03:04
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-28 12:18:28
 */

#include "lineage.h"

ostream &operator<<(ostream &os, const Lineage &lng) {
  vector<string> ranks;
  separateLineage(lng.name, ranks);
  for(auto& str : ranks){
    str = lastNameNoRank(str);
  }
  os << "[\"";
  os << ranks.back() << "\",\"";
  os << strjoin(ranks.begin(), ranks.end()-1, "\",\"");
  os << "\"]";
  return os;
};

/********************************************************************************
 * @brief Construct a new Lineage Handle:: Lineage Handle object.
 *        basic option of Lineage
 * @param taxfile
 ********************************************************************************/
LineageHandle::LineageHandle() { rank = rank->create(); };

LineageHandle::LineageHandle(const string &dbpath, const string &taxfile,
                             const string &revfile)
    : tdbpath(dbpath), taxfile(taxfile), revfile(revfile) {
  rank = rank->create();
};

/********************************************************************************
 * @brief obtained the lineage string from the taxfile
 *
 * @param lngs a list of lineage for query and retrun
 ********************************************************************************/
void LineageHandle::getLineage(vector<Lineage> &lngs) {

  // find the lineage in the taxfile
  TaxMap taxmap;
  readTaxFile(taxmap);
  size_t nHit(0);
  if (!taxmap.empty()) {
    for (auto &nm : lngs) {
      auto iter = taxmap.find(nm.name);
      if (iter != taxmap.end()) {
        nm.name = rank->lineage(iter->second, nm.name);
        nHit++;
        nm.def = true;
      }
    }
    theInfo("Found " + to_string(nHit) + "/" + to_string(lngs.size()) +
            " lineages in file " + taxfile);
  }

  // find lineage in database
  if (nHit < lngs.size() && fileExists(tdbpath)) {
    TaxaDB taxadb(tdbpath);
    size_t mHit(0);
    for (auto &lng : lngs) {
      if (!lng.def) {
        string res = taxadb.search(lng.name);
        if (!res.empty()) {
          lng.name = rank->lineage(res, lng.name);
          mHit++;
        }
      }
    }
    theInfo("Found " + to_string(mHit) + "/" + to_string(lngs.size()) +
            " lineages in taxa database " + tdbpath);
  }

  // format the lineage and do revsion
  for (auto &lng : lngs) {
    rank->format(lng.name);
  }

  // do revision
  Revision rev(revfile);
  if (!rev.empty()) {
    for (auto &lng : lngs) {
      rev.revise(lng.name);
    }
  }

  // check Repeat names
  checkRepeats(lngs);

  // set the unclassified items
  for (auto &nm : lngs) {
    nm.def = rank->wellDefined(nm.name);
  }
};

// read lineage string from input taxon file
void LineageHandle::readTaxFile(TaxMap &taxmap) {

  // get lineage map from file
  ifstream tax(taxfile);
  if (!tax) {
    theInfo("Open file " + taxfile +
            " failed, all lineages will obtained from taxon database");
  } else {
    // read the taxonomy map
    regex matchReg("([^; ]+)[; ]+([^; ]+)");
    for (string line; getline(tax, line);) {
      if (!line.empty()) {
        smatch matchs;
        regex_search(line, matchs, matchReg);
        if (!matchs.empty()) {
          taxmap[matchs[1].str()] = matchs[2].str();
        } else {
          theInfo("Error parse in " + line + " in file " + taxfile);
        }
      }
    }
  }
  tax.close();
};

// check the repeat items
void LineageHandle::checkRepeats(vector<Lineage> &lngs) {
  map<string, int> nameMap;
  for (auto &nm : lngs) {
    if (nameMap.find(nm.name) == nameMap.end()) {
      nameMap[nm.name] = 1;
    } else {
      theInfo("Found repeat strain: " + nm.name);
      nm.name += ".repeat-" + to_string(++nameMap[nm.name]);
    }
  }
}
