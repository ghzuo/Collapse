/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-05-07 08:57:21
 */

#include "collapse.h"

void collapse(int argc, char *argv[]) {

  // get the input arguments
  RunArgs myargs(argc, argv);

  /************************************************************************
   ******* read the tree and checked *************************************/
  Node *aTree = new Node;
  aTree->innwk(myargs.infile);

  /**********************************************************************
   ********* set for the lineage system and get lineage *****************/
  TaxaRank *rank = TaxaRank::create();
  rank->initial(myargs.rankfile, myargs.outrank);
  LngData lngs(myargs.taxadb, myargs.taxfile, myargs.taxrev);

  // get the leafs lineage
  vector<Node *> allLeafs;
  aTree->getLeafs(allLeafs);
  theInfo("There are " + to_string(allLeafs.size()) +
          " leafs in the phylogenetic tree");

  // set lineage of leafs
  vector<string> nmlist;
  for (auto nd : allLeafs) {
    nmlist.emplace_back(nd->name);
  }
  lngs.getLineage(nmlist);

  /*************************************************************************
   *********** set lineage of Nodes and rooting unrooted tree **************/
  // set lineage of leaf
  for (int i = 0; i < allLeafs.size(); ++i) {
    allLeafs[i]->setOneLeaf(lngs[i].name, lngs[i].def);
  }

  // annotate the branches
  if (aTree->children.size() == 2) {
    // update the rooted tree
    aTree->updateRootedTree();
  } else if (!myargs.outgrp.empty()) {
    // root the unrooted tree by input outgroup name
    Node *result = aTree->rootingByOutgrp(myargs.outgrp);
    if (result == NULL) {
      theInfo("Rooting the tree by taxonomy");
      aTree = aTree->rootingByTaxa();
    } else {
      aTree = result;
    }
  } else {
    // rooting tree by taxonomy
    aTree = aTree->rootingByTaxa();
    aTree->balanceTree(myargs.rootMeth, myargs.otuLevel);
  }

  theInfo("Annotated all nodes of tree by lineages");

  /**************************************************************************
   ************ get the statistic of Taxonomy and annotate tree *************/
  // get the statics of the two taxa system (defined and undefined)
  Taxa aTaxa(lngs);

  // annotate taxonomy of nodes
  aTaxa.annotate(aTree);

  theInfo("Done statistics of the taxonomy");

  /***************************************************************************
   *********  output data ****************************************************/
  if (myargs.forApp) {
    out4app(lngs, aTaxa, aTree, myargs.outPref + ".json");
  } else if (myargs.forWeb) {
    out4serv(lngs, aTaxa, aTree, myargs);
  } else {
    output(lngs, aTaxa, aTree, myargs);
  }
}
/****************************************************************************
 ******************************  End main program ***************************
 ****************************************************************************/

RunArgs::RunArgs(int argc, char **argv)
    : infile(""), taxrev(""), outgrp(""), taxfile(""), forWeb(false),
      forApp(false), predict(false), itol(false), lngfile(""), clevel(""),
      rootMeth("mad"), otuLevel("") {

  program = argv[0];
  string outname("collapsed");
  string supdir("./");
  string kstr;

  char ch;
  while ((ch = getopt(argc, argv, "i:d:D:o:m:r:t:s:T:R:O:L:l:C:a:u:IWPAJqh")) !=
         -1) {
    switch (ch) {
    case 'i':
      infile = optarg;
      break;
    case 'D':
      supdir = optarg;
      break;
    case 'd':
      taxadb = optarg;
      break;
    case 'o':
      outname = optarg;
      break;
    case 'm':
      taxrev = optarg;
      break;
    case 'r':
      outrank = optarg;
      break;
    case 'R':
      rankfile = optarg;
      break;
    case 'O':
      outgrp = optarg;
      break;
    case 'l':
      taxfile = optarg;
      break;
    case 'L':
      lngfile = optarg;
      break;
    case 'C':
      clevel = optarg;
      break;
    case 'a':
      rootMeth = toLower(optarg);
      break;
    case 'u':
      otuLevel = toUpper(optarg);
      break;
    case 'I':
      itol = true;
      break;
    case 'W':
      forWeb = true;
      break;
    case 'A':
      forApp = true;
      break;
    case 'P':
      predict = true;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'J':
      theJson.format = true;
      break;
    case 'h':
      usage();
    case '?':
      usage();
    }
  }

  // the workdir
  addsuffix(supdir, '/');

  // the input tree file
  if (infile.empty())
    infile = supdir + "Tree.nwk";

  // the input lineage file
  if (taxfile.empty()) {
    taxfile = supdir + "Lineage.lsn";
    if (!fileExists(taxfile))
      taxfile = supdir + "Lineage.csv";
  }

  // the output prefix
  outPref = supdir + outname;

  // for the default
  if (taxadb.empty()) {
    taxadb = supdir + "taxadb.gz";
    if (!fileExists(taxadb)) {
      taxadb = supdir + "taxdump.tar.gz";
      if (!fileExists(taxadb)) {
        taxadb = supdir + "taxdump/";
      }
    }
  }
}

void RunArgs::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -D ./ ]              The work directory, default: ./\n"
       << " [ -i Tree.nwk ]        Input newick tree, default: Tree.nwk\n"
       << " [ -o collapsed ]       Set prefix name of output files, \n"
       << "                        default: collapsed\n"
       << " [ -m <Revision.txt> ]  Lineage revision file for batch edit,\n"
       << "                        default: None\n"
       << " [ -l Lineage.txt ]     Input lineage file for leaves of tree, \n"
       << "                        default: Lineage.txt or Lineage.csv\n"
       << " [ -d taxadb.gz ]       Taxonomy data file or directory,\n"
       << "                        default: taxadb.gz or taxdump.tar.gz\n"
       << " [ -R <None> ]          List of rank names and abbravitions,\n"
       << "                        default: set by program\n"
       << " [ -r DKPCOFGS ]        Abbreviations of output taxon rank,\n"
       << "                        default: set by program\n"
       << " [ -C <None> ]          Collapse tree on the taxon level,\n"
       << "                        default: the top division taxon level\n"
       << " [ -O <Outgroup> ]      Set the outgroup for the unroot tree.\n"
       << "                        default: None, rearranged by taxonomy\n"
       << " [ -I ]                 Output newick and annotate files for itol\n"
       << "                        default: No\n"
       << " [ -P ]                 Output prediction for undefined leafs\n"
       << " [ -q ]                 Run command in quiet mode\n"
       << " [ -h ]                 Display this information\n"
       << endl;
  exit(1);
}

/**************************************************************************
 * @brief output the tree (nwk and/or json), and statistics of taxonomy
 *
 * @param aTaxa the statisics of taxonomy
 * @param aTree the tree
 * @param myargs the output file name and switches
 ***************************************************************************/
void output(const LngData &lngs, Taxa &aTaxa, Node *aTree, RunArgs &myargs) {

  /// output the monophyly
  aTaxa.outStatitics(myargs.outPref + ".unit");

  /// output the entropy
  aTaxa.outEntropy(myargs.outPref + ".entropy");

  /// output the annotated newick file
  if (myargs.itol) {
    aTree->outitol(myargs.outPref + "-iToL_tree.nwk");
    vector<Node *> nodes;
    aTree->getDescendants(nodes);
    outItolNodes(nodes, myargs.outPref + "-iToL_nodelist.txt");
    outItolLabel(nodes, myargs.outPref + "-iToL_label.txt");
    outItolPopup(aTaxa, nodes, myargs.outPref + "-iToL_popup.txt");
    outItolSymbol(aTaxa, nodes, myargs.outPref + "-iToL_symbol.txt");

    // get division rank
    set<string> division;
    getDivision(aTaxa, myargs.clevel, division);
    if (division.size() > 1) {
      theInfo("The phylogenetic tree will shown in " +
              to_string(division.size()) + " classes");
      outItolStrap(nodes, division, myargs.outPref + "-iToL_strap.txt");
      outItolCollapse(nodes, division, myargs.outPref + "-iToL_collapse.txt");
    }
  } else {
    aTree->outnwk(myargs.outPref + ".nwk");
  }

  // for undefined items
  if (myargs.predict && aTree->nxleaf > 0) {
    /// output the unclassified items
    // vector<string> strName;
    // aTree->getUndefineNames(strName);
    // aTaxa.outUnclass(strName, myargs.outPref + ".unclass");

    /// output the prediction of unclassifed items
    aTree->outPrediction(myargs.outPref + ".predict");
  }

  // oupt lineage of leafs
  ofstream leaf(myargs.outPref + ".lineage");
  lngs.outcsv(leaf);
  leaf.close();
};

// output for web sever
void out4serv(const LngData &lngs, Taxa &aTaxa, Node *aTree, RunArgs &myargs) {
  /// open the ostream for tree page
  stringstream jstree;
  outTreeJson(aTaxa, aTree, jstree);
  theJson(myargs.outPref + "-phylogeny.json", jstree.str());

  /// open the ostream for taxa page
  stringstream jstaxa;
  outTaxaJson(aTaxa, aTree, jstaxa);
  theJson(myargs.outPref + "-taxonomy.json", jstaxa.str());

  /// output put lineage file
  if (!myargs.lngfile.empty()) {
    lngs.output(myargs.lngfile, "json");
  }
};

// output the complex json for cltree app
void out4app(const LngData &lngs, Taxa &aTaxa, Node *aTree,
             const string &file) {
  /// open the ostream
  stringstream jsonbuf;

  // output the phylogeny
  jsonbuf << "{\"phylogeny\":";
  outTreeJson(aTaxa, aTree, jsonbuf);
  jsonbuf << ", " << endl;

  // output the taxonomy
  jsonbuf << "\"taxonomy\":";
  outTaxaJson(aTaxa, aTree, jsonbuf);
  jsonbuf << "," << endl;

  // output the lineage
  jsonbuf << "\"lineage\":";
  lngs.outjson(jsonbuf);
  jsonbuf << "}" << endl;

  /// output the taxonomy list with relationship
  theJson(file, jsonbuf.str());
}

void outTaxaJson(Taxa &aTaxa, Node *aTree, ostream &os) {
  // output the taxon system
  os << "{\"taxa\":";
  if (aTree->nleaf == 0) {
    os << "[]";
  } else {
    aTaxa.outJsonTax(os);
  }
  os << "," << endl;

  // output unclassified info
  if (aTree->nxleaf > 0) {
    /// output the unclassified items
    vector<string> strName;
    aTree->getUndefineNames(strName);
    os << "\"unclass\":";
    aTaxa.outJsonUnclass(strName, os);
    os << "," << endl;
  }
  // output the entropy
  os << "\"rank\":";
  aTaxa.rank->outRanksJson(os);
  os << "}";
};

void outTreeJson(Taxa &aTaxa, Node *aTree, ostream &os) {
  //... check the upload genomes
  aTree->checkUploaded();
  os << "{\"tree\":";
  aTree->outjson(os);
  os << "," << endl;

  // output the level statistics
  os << "\"statistics\":";
  aTaxa.outJsonEntropy(os);
  os << "}";
};

/********************************************************************************
 * @brief the function for output iToL files
 *
 * @param os
 * @param file
 * @param type
 ********************************************************************************/
void outItolNodes(const vector<Node *> &nodes, const string &file) {
  ofstream os(file);
  for (auto &nd : nodes) {

    if (nd->isLeaf()) {
      os << nd->id << "\t"
         << "leaf"
         << "\t" << nd->name.substr(nd->name.find_first_of('|') + 1);
    } else {
      os << "I" << nd->id << "\t";
      int nlvl = nd->taxLevel - nd->parent->taxLevel;
      if (nd->nleaf > 0 && nlvl > 0) {
        string lngstr = nd->name;
        lngstr.erase(remove(lngstr.begin(), lngstr.end(), '|'), lngstr.end());
        vector<string> lngvec;
        separateLineage(lngstr, lngvec);
        os << "clade"
           << "\t";
        for (int i = lngvec.size() - nlvl; i < lngvec.size(); ++i)
          os << lngvec[i];
      } else {
        os << "-"
           << "\t" << nd->name.substr(nd->name.find_first_of('|') + 1);
      }
    }

    os << endl;
  }
};

void outItolLabel(const vector<Node *> &nodes, const string &file) {
  ofstream os(file);
  itolHeader(os, file, "LABELS");

  // the data section
  os << "DATA\n" << endl;
  for (auto &nd : nodes) {
    if (nd->isLeaf()) {
      os << nd->id << "\t" << lastNameNoRank(nd->name) << endl;
    } else if (nd->taxLevel > nd->parent->taxLevel) {
      os << "I" + to_string(nd->id) << "\t" << lastNameNoRank(nd->name) << endl;
    }
  }

  os.close();
};

void outItolPopup(const Taxa &aTaxa, const vector<Node *> &nodes,
                  const string &file) {
  ofstream os(file);
  itolHeader(os, file, "POPUP_INFO");

  // the data section
  os << "DATA\n" << endl;
  for (auto &nd : nodes) {
    if (nd->isLeaf()) {
      os << nd->id << "\tSpecies information:\t" << itolPopusStr(nd, aTaxa)
         << endl;
    } else if (nd->nleaf > 0 && nd->taxLevel > nd->parent->taxLevel) {
      os << "I" + to_string(nd->id) << "\tClade information:\t"
         << itolPopusStr(nd, aTaxa) << endl;
    }
  }

  os.close();
};

void outItolStrap(const vector<Node *> &nodes, const set<string> &division,
                  const string &file) {
  // the file header
  ofstream os(file);
  itolHeader(os, file, "DATASET_COLORSTRIP");
  os << "DATASET_LABEL\tcolor_strip\n"
     << "COLOR_BRANCHES\t1\n"
     << "STRIP_WIDTH\t25\n"
     << "MARGIN\t0\n"
     << "BORDER_WIDTH\t0\n"
     << "SHOW_INTERNAL\t0\n";

  // the data section
  os << "DATA\n" << endl;

  // the color map
  int theRank = nRanks(*division.begin()) - 1;
  map<string, string> colorMap;
  for (auto &d : division) {
    colorMap[d] = "X";
  }
  getColorMap(colorMap);

  for (auto &nd : nodes) {
    string lngstr = nd->name;
    lngstr.erase(remove(lngstr.begin(), lngstr.end(), '|'), lngstr.end());
    vector<string> nmlist;
    parseLineage(lngstr, nmlist);
    if (nmlist.size() > theRank) {
      auto iter = colorMap.find(nmlist[theRank]);
      if (iter != colorMap.end()) {
        if (!nd->isLeaf())
          os << "I";
        os << to_string(nd->id) << "\t" << iter->second << "\t"
           << lastNameNoRank(iter->first) << endl;
      }
    }
  }

  os.close();
};

void outItolCollapse(const vector<Node *> &nodes, const set<string> &division,
                     const string &file) {
  // the file header
  ofstream os(file);
  itolHeader(os, file, "COLLAPSE");

  // the data section
  os << "DATA\n" << endl;

  int theRank = nRanks(*division.begin());
  for (auto &nd : nodes) {
    if (!nd->isLeaf()) {
      if (nd->parent->taxLevel < theRank && nd->taxLevel >= theRank) {
        os << "I" << nd->id << endl;
      }
    }
  }
  os.close();
};

void outItolSymbol(const Taxa &aTaxa, const vector<Node *> &nodes,
                   const string &file) {
  // the file header
  ofstream os(file);
  itolHeader(os, file, "DATASET_SYMBOL");
  os << "DATASET_LABEL\tsymbol\n"
     << "MAXIMUM_SIZE\t20\n"
     << "COLOR\t#ffff00\n";

  // the data section
  os << "DATA\n" << endl;

  for (auto &nd : nodes) {
    if (!nd->isLeaf() && nd->nleaf > 0 && nd->taxLevel > nd->parent->taxLevel) {
      os << "I" << to_string(nd->id) << "\t" << 3 << "\t" << 30 << "\t"
         << "#999999"
         << "\t" << 1 << "\t" << 1.0 << endl;
    }
  }
  os.close();
};

void itolHeader(ostream &os, const string &file, const string &type) {
  // test output file stream
  if (!os) {
    cerr << "Open " << file << " for write failed" << endl;
    exit(3);
  }

  // the title
  os << type << "\n"
     << "SEPARATOR TAB\n";
};

void getDivision(const Taxa &aTaxa, const string &str, set<string> &division) {
  int theRank = 0;
  if (!str.empty()) {
    theRank = aTaxa.rank->rankindex(str);
    if (theRank == 0) {
      theInfo("Cann't find the rank level: " + str);
    }
  }

  if (theRank != 0) {
    for (auto &tax : aTaxa.def.state) {
      if (nRanks(tax.first) == theRank)
        division.insert(tax.first);
    }
  } else {
    theInfo("Use the top division rank level");
    getTopDivision(aTaxa, division);
  }
};

void getTopDivision(const Taxa &aTaxa, set<string> &division) {
  vector<set<string>> rankset(aTaxa.rank->outrank.size() + 1);
  for (auto &tax : aTaxa.def.state) {
    rankset[nRanks(tax.first)].insert(tax.first);
  }

  for (auto &rk : rankset) {
    if (rk.size() > 1) {
      division = rk;
      break;
    }
  }
}

void getColorMap(map<string, string> &colorMap) {
  float delta = 360 / colorMap.size();
  int i = 0;
  for (auto &cm : colorMap) {
    vector<int> cv{(int)(delta * i++), 99, 90};
    hsv2rgb(cv);

    stringstream buf;
    buf << "rgba(";
    for (auto c : cv)
      buf << c << ",";
    buf << "1.0)";
    cm.second = buf.str();
  }
}

string itolPopusStr(Node *nd, const Taxa &aTaxa) {

  stringstream buf;

  buf << "<div class='tPop'>"
      << "<h1>" << lastNameNoRank(nd->name) << "</h1>"
      << "<h2>Branch length: " << nd->length << "</h2>"
      << "<h2>" << nd->nleaf;
  if (nd->nxleaf > 0)
    buf << "+" << nd->nxleaf;
  buf << " Leaves</h2>"
      << "<br><h1>Lineage Information: </h1>";

  string lngstr = nd->name;
  lngstr.erase(remove(lngstr.begin(), lngstr.end(), '|'), lngstr.end());
  vector<string> nmlist;
  parseLineage(lngstr, nmlist);
  size_t nOutput = nmlist.size() < aTaxa.rank->outrank.size()
                       ? nmlist.size()
                       : aTaxa.rank->outrank.size();

  buf << "<table class = 'lineage'>";
  for (size_t i = 0; i < nOutput; ++i) {
    buf << "<tr>"
        << "<th>" << aTaxa.rank->outrank[i].first << ": </th>"
        << "<td>" << lastNameNoRank(nmlist[i]);

    auto iter = aTaxa.def.state.find(nmlist[i]);
    buf << "{" << iter->second.nStrain << ", " << iter->second.distract.size()
        << "}</td></tr>";
  }
  buf << "</table>";

  buf << "</div>";
  return buf.str();
};