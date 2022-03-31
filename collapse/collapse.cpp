/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-31 19:29:17
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
    if(myargs.withNHX)
      aTree->setNwkWithNHX();
    output(lngs, aTaxa, aTree, myargs);
  }
}
/****************************************************************************
 ******************************  End main program ***************************
 ****************************************************************************/

RunArgs::RunArgs(int argc, char **argv)
    : infile(""), taxrev(""), outgrp(""), taxfile(""), forWeb(false),
      forApp(false), predict(false), withNHX(false), lngfile("") {

  program = argv[0];
  string outname("collapsed");
  string supdir("./");
  string kstr;

  char ch;
  while ((ch = getopt(argc, argv, "i:d:D:o:m:t:s:T:R:O:L:l:NWPAJqh")) != -1) {
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
    case 'N':
      withNHX = true;
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
       << " [ -O <Outgroup> ]      Set the outgroup for the unroot tree.\n"
       << "                        default: None, rearranged by taxonomy\n"
       << " [ -N ]                 Output newick with NHX for display it\n"
       << "                        in iToL Website, default: No\n"
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
  aTree->outnwk(myargs.outPref + ".nwk");

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
