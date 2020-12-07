/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-01-02 15:42:07
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-03 14:03:32
 */

#include "taxsys.h"

/********************************************************************************
 * default value for static members of taxsys
 ********************************************************************************/
const string TaxSys::rootTaxon("<B>Cellular_Organisms");
Abbr TaxSys::division{{"Domain", 'D'}, {"Kingdom", 'K'}, {"Phylum", 'P'},
                      {"Class", 'C'},  {"Order", 'O'},   {"Family", 'F'},
                      {"Genus", 'G'},  {"Species", 'S'}};
void TaxSys::setDivision(const Abbr &dv) { division = dv; }

/********************************************************************************
 * @brief Construct a new Tax Sys:: Tax Sys object
 *
 * @param names the lineage of a type leafs (defined or undefined)
 ********************************************************************************/
TaxSys::TaxSys(const vector<string> &names) { initial(names); };

void TaxSys::initial(const vector<string> &names) {
  //set the total number of the taxonomy set
  nStrain = names.size();

  // set the state map
  for (auto &nm : names) {
    vector<string> taxlist;
    parseLineage(nm, taxlist);
    taxlist.pop_back();
    for (auto &tax : taxlist) {
      if (state.find(tax) != state.end()) {
        ++state[tax].nStrain;
      } else {
        state[tax] = TaxonState();
      }
    }
  }
};

void TaxSys::annotateBranch(size_t nleaf, int nClade, string &name,
                            size_t &taxSize) {

  if (name.empty()) {
    name = "|" + rootTaxon;
    taxSize = nStrain;
  } else {
    // reset the name of the node
    vector<string> taxlist;
    parseLineage(name, taxlist);

    // set the clades
    for (size_t i = 1; i <= nClade; ++i) {
      state[taxlist[taxlist.size() - i]].distract.emplace_back(nleaf);
    }

    // set taxSize
    auto iter = taxlist.rbegin();
    auto stateIter = state.find(*iter);
    taxSize = stateIter->second.nStrain;
    if (stateIter->second.nStrain == nleaf) {
      stateIter->second.monophy = true;
    }

    // set the name of node
    size_t npos(0);
    for (++iter; iter != taxlist.rend(); ++iter) {
      stateIter = state.find(*iter);
      if (stateIter->second.nStrain == nleaf) {
        stateIter->second.monophy = true;
      } else {
        npos = iter->size();
        break;
      }
    }
    name.insert(npos, "|");
  }
};

void TaxSys::annotateLeaf(int nClade, string &name, size_t &taxSize) {

  // for the case only one strain
  taxSize = 1;

  // reset the name of the node
  vector<string> taxlist;
  parseLineage(name, taxlist);
  taxlist.pop_back();

  // set the clades
  for (int i = 1; i < nClade; ++i) {
    state[taxlist[taxlist.size() - i]].distract.emplace_back(1);
  }

  // set node name
  size_t npos(0);
  for (auto iter = taxlist.rbegin(); iter != taxlist.rend(); ++iter) {
    auto stateIter = state.find(*iter);
    if (stateIter->second.nStrain == 1) {
      stateIter->second.monophy = true;
    } else {
      npos = iter->size();
      break;
    }
  }
  name.insert(npos, "|");
};

//..... Combine the classified and unclassified tax system
/********************************************************************************
 * @brief Construct a new Taxa:: Taxa object
 *
 * @param nmlist a list of the lineages for all leafs
 ********************************************************************************/
Taxa::Taxa(const vector<Lineage> &lngs) {
  vector<string> defNames;
  vector<string> undefNames;
  for (auto &lng : lngs) {
    if (lng.def) {
      defNames.emplace_back(lng.name);
    } else {
      undefNames.emplace_back(lng.name);
    }
  }

  def.initial(defNames);
  undef.initial(undefNames);
};

void Taxa::annotate(Node *aTree) {
  // set taxonomy of leafs
  vector<Node *> allLeafs;
  aTree->getLeafs(allLeafs);
  for (auto &nd : allLeafs) {
    if (nd->unclassified) {
      undef.annotateLeaf(nd->nClade(), nd->name, nd->taxSize);
    } else {
      def.annotateLeaf(nd->nClade(), nd->name, nd->taxSize);
    }
  }

  // annotetate taxonomy of all nodes
  vector<Node *> allBraches;
  aTree->getBranches(allBraches);
  for (auto &nd : allBraches) {
    if (nd->unclassified) {
      undef.annotateBranch(nd->nxleaf, nd->nClade(), nd->name, nd->taxSize);
    } else {
      def.annotateBranch(nd->nleaf, nd->nClade(), nd->name, nd->taxSize);
    }
  }
};

/********************************************************************************
 * @brief output the list according the lineage for web server
 *
 * @param os
 ********************************************************************************/
void Taxa::outTax(const string &file) {
  ofstream os(file);
  if (!os) {
    cerr << "Open " << file << " for write failed" << endl;
    exit(3);
  }

  outTax(os);
  os.close();
};

void Taxa::outTax(ostream &os) {
  for (const auto &st : def.state) {
    os << lastName(st.first) << ":" << st.second.nStrain;
    if (st.second.monophy)
      os << " ++" << endl;
    else
      os << " --" << endl;
  }
};

/********************************************************************************
 * @brief output the unclassified items
 *
 * @param strName
 * @param os
 ********************************************************************************/
void Taxa::outUnclass(vector<string> &strName, const string &file) {
  ofstream os(file);
  if (!os) {
    cerr << "Open " << file << " for write failed" << endl;
    exit(3);
  }

  outUnclass(strName, os);
  os.close();
};

void Taxa::outUnclass(vector<string> &strName, ostream &os) {
  vector<string> taxName;
  for (auto &str : undef.state)
    taxName.push_back(str.first);

  for (auto &str : strName)
    str.erase(str.find('|'), 1);
  sort(strName.begin(), strName.end());

  vector<string> names;
  merge(taxName.begin(), taxName.end(), strName.begin(), strName.end(),
        back_inserter(names));
  for (auto &str : names)
    os << lastName(str) << ':' << undef.state[str].nStrain << endl;
};

/********************************************************************************
 * @brief output the entropy
 *
 * @param os
 ********************************************************************************/
void Taxa::outEntropy(const string &file) {
  ofstream os(file);
  if (!os) {
    cerr << "Open " << file << " for write failed!" << endl;
    exit(3);
  }

  outEntropy(os);
  os.close();
};

void Taxa::outEntropy(ostream &os) {

  map<char, TaxLevelState> taxlev;
  for (const auto &atax : TaxSys::division) {
    TaxLevelState tl;
    taxlev[atax.second] = tl;
  }

  for (auto &stateIter : def.state) {

    string theItem = lastName(stateIter.first);
    char abbrT = theItem[1];

    // get the collapse state and sigle strain taxon
    if (stateIter.second.nStrain == 1) {
      ++taxlev[abbrT].nSolo;
    } else {
      // get number of non-solo taxon
      if (stateIter.second.monophy)
        ++taxlev[abbrT].nMono;
      else
        ++taxlev[abbrT].nPoly;

      // get the entropy
      taxlev[abbrT].sTax +=
          stateIter.second.nStrain * log2(double(stateIter.second.nStrain));
      for (auto &n : stateIter.second.distract) {
        taxlev[abbrT].sTree += n * log2(double(n));
      }
    }
  }

  double maxEntropy = log2(double(def.nStrain));
  for (const auto &atax : TaxSys::division) {
    int nTotal = taxlev[atax.second].nSolo + taxlev[atax.second].nMono +
                 taxlev[atax.second].nPoly;
    double sTax = maxEntropy - taxlev[atax.second].sTax / def.nStrain;
    double sTree = maxEntropy - taxlev[atax.second].sTree / def.nStrain;

    // output taxlevel data
    os << atax.first << "\t" << taxlev[atax.second].nSolo << "\t"
       << taxlev[atax.second].nMono << "\t" << nTotal << "\t" << sTax << "\t"
       << sTree << endl;
  }

  // output the strain data
  os << "Strain"
     << "\t" << def.nStrain << "\t" << 0 << "\t" << def.nStrain << "\t"
     << maxEntropy << "\t" << maxEntropy << endl;
};

/********************************************************************************
 * @brief output the statistics according taxonomy levels
 *
 * @param os
 ********************************************************************************/
void Taxa::outStatitics(const string &file) {
  ofstream os(file);
  if (!os) {
    cerr << "Open " << file << " for write failed" << endl;
    exit(3);
  }

  outStatitics(os);
  os.close();
};

void Taxa::outStatitics(ostream &os) {
  for (const pair<string, char> &atax : TaxSys::division) {

    // get the abbr
    char abbr = atax.second;

    // get the statistic of monophyly
    size_t nMul(0), nSol(0);
    vector<string> tlist;
    for (auto &stateIter : def.state) {
      string theItem = lastName(stateIter.first);
      if (theItem[1] == abbr) {
        // get the number of strain
        theItem.append(":");
        theItem.append(to_string(stateIter.second.nStrain));

        // get the collapse state and sigle strain taxon
        if (stateIter.second.monophy) {
          theItem.append("\t++");
          if (stateIter.second.nStrain == 1)
            ++nSol;
          else
            ++nMul;
        } else
          theItem.append("\t--");

        // get the divide
        theItem.append("\t{");
        theItem.append(to_string(stateIter.second.distract[0]));

        for (int i = 1; i < stateIter.second.distract.size(); ++i) {
          theItem.append("|");
          theItem.append(to_string(stateIter.second.distract[i]));
        }
        theItem.append("}");
        tlist.emplace_back(theItem);
      }
    }

    // output the result
    os << atax.first << " (" << nSol << "+" << nMul << "/" << tlist.size()
       << "):" << endl;

    sort(tlist.begin(), tlist.end());
    for (auto &str : tlist)
      os << str << endl;
    os << endl;
  }
};
