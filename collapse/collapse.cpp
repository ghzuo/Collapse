/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-09-01 12:54:33
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2021-01-18 14:17:48
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
  rank->initial(myargs.taxfile, myargs.abfile, myargs.abtype);
  LineageHandle lngHandle(myargs.taxadb, myargs.taxfile, myargs.taxrev);

  // get the leafs lineage
  vector<Node *> allLeafs;
  aTree->getLeafs(allLeafs);
  vector<Lineage> lngs;
  for (auto nd : allLeafs) {
    lngs.emplace_back(nd->name);
  }
  theInfo("There are " + to_string(allLeafs.size()) +
          " leafs in the phylogenetic tree");

  // set lineage of leafs
  lngHandle.getLineage(lngs);

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
    output(lngs, aTaxa, aTree, myargs);
  }
}
/****************************************************************************
 ******************************  End main program ***************************
 ****************************************************************************/

RunArgs::RunArgs(int argc, char **argv)
    : infile(""), taxrev(""), outgrp(""), taxfile(""), forWeb(false),
      forApp(false), predict(true) {

  program = argv[0];
  string outname("collapsed");
  string supdir("./");
  string kstr;

  char ch;
  while ((ch = getopt(argc, argv, "i:d:D:o:r:t:s:T:R:O:WPAqh")) != -1) {
    switch (ch) {
    case 'i':
      infile = optarg;
      break;
    case 'd':
      supdir = optarg;
      break;
    case 'D':
      taxadb = optarg;
      break;
    case 'o':
      outname = optarg;
      break;
    case 'r':
      taxrev = optarg;
      break;
    case 't':
      abtype = optarg;
      break;
    case 'R':
      abfile = optarg;
      break;
    case 'O':
      outgrp = optarg;
      break;
    case 'T':
      taxfile = optarg;
      break;
    case 'W':
      forWeb = true;
      break;
    case 'A':
      forApp = true;
      break;
    case 'P':
      predict = false;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'h':
      usage();
    case '?':
      usage();
    }
  }

  // the workdir
  addsuffix(supdir, '/');

  // the input items
  if (infile.empty())
    infile = supdir + "Tree.nwk";
  if (taxrev.empty())
    taxrev = supdir + "Lineage.rev";
  if (taxfile.empty())
    taxfile = supdir + "Lineage.list";

  // the output prefix
  outPref = supdir + outname;

  // for the default
  if (taxadb.empty()) {
    taxadb = supdir + "taxadb.gz";
    if (!fileExists(taxadb)) {
      taxadb = supdir + "taxdump.tar.gz";
      if (!fileExists(taxadb)) {
        taxadb = supdir + "taxdump";
      }
    }
  }
}

void RunArgs::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -d ./ ]            The work directory, default: ./\n"
       << " [ -i Tree.nwk ]      Input newick tree file, default: Tree.nwk\n"
       << " [ -o collapsed ]     Output prefix name: default: collapsed\n"
       << " [ -r Lineage.rev ]   Lineage revise file for batch edit,\n"
       << "                      default: Lineage.rev\n"
       << " [ -T Lineage.list ]  Lineage file for leafs of tree, \n"
       << "                      default: Lineage.list\n"
       << " [ -D taxadb.gz ]     Taxa database file or directory,\n"
       << "                      default: taxadb.gz or taxdump/\n"
       << " [ -R <None> ]        Abbravition list for taxon rank name,\n"
       << "                      default: by program\n"
       << " [ -t DKPCOFGS ]      Abbreviation of output taxon rank,\n"
       << "                      default: DKPCOFGS\n"
       << " [ -O <Outgroup> ]    Set the outgroup for the unroot tree.\n"
       << "                      default: None, rearranged by taxonomy\n"
       << " [ -P ]               Output the prediction for undefined leafs\n"
       << " [ -q ]               Run command in quiet mode\n"
       << " [ -h ]               Display this information\n"
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
void output(const vector<Lineage> &lngs, Taxa &aTaxa, Node *aTree,
            RunArgs &myargs) {

  /// output the monophyly
  aTaxa.outStatitics(myargs.outPref + ".unit");

  /// output the entropy
  aTaxa.outEntropy(myargs.outPref + ".entropy");

  /// output the annotated newick file
  aTree->outnwk(myargs.outPref + "-annotated.nwk");

  // for undefined items
  if (aTree->nxleaf > 0) {
    /// output the unclassified items
    vector<string> strName;
    aTree->getUndefineNames(strName);
    aTaxa.outUnclass(strName, myargs.outPref + ".unclass");

    /// output the prediction of unclassifed items
    if (myargs.predict) {
      aTree->outPrediction(myargs.outPref + ".predict");
    }
  }

  // oupt lineage of leafs
  ofstream leaf(myargs.outPref + ".lineage");
  for (auto &lng : lngs) {
    leaf << lng.name << endl;
  }
  leaf.close();
};

// output for web sever
void out4serv(const vector<Lineage> &lngs, Taxa &aTaxa, Node *aTree,
              RunArgs &myargs) {
  /// open the ostream for tree page
  ofstream TREE(myargs.outPref + "-phylogeny.json");
  outTreeJson(aTaxa, aTree, TREE);
  TREE.close();

  /// open the ostream for taxa page
  ofstream TAXA(myargs.outPref + "-taxonomy.json");
  outTaxaJson(aTaxa, aTree, TAXA);
  TAXA.close();

  /// open the ostream for lineage of leafs
  ofstream LNGS(myargs.outPref + "-lineage.json");
  outLngsJson(lngs, aTaxa, LNGS);
  LNGS.close();
};

// output the complex json for cltree app
void out4app(const vector<Lineage> &lngs, Taxa &aTaxa, Node *aTree,
             const string &file) {
  /// open the ostream
  ofstream JSON(file);

  // output the phylogeny
  JSON << "{\"phylogeny\":";
  outTreeJson(aTaxa, aTree, JSON);
  JSON << ", " << endl;

  // output the taxonomy
  JSON << "\"taxonomy\":";
  outTaxaJson(aTaxa, aTree, JSON);
  JSON << "," << endl;

  // output the lineage
  JSON << "\"lineage\":";
  outLngsJson(lngs, aTaxa, JSON);
  JSON << "}" << endl;

  /// output the taxonomy list with relationship
  JSON.close();
}

void outTaxaJson(Taxa &aTaxa, Node *aTree, ostream &os) {
  // output the taxon system
  os << "{\"taxa\":";
  aTaxa.outJsonTax(os);
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
  os << "}" << endl;
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
  os << "}" << endl;
};

void outLngsJson(const vector<Lineage> &lngs, Taxa &aTaxa, ostream &os) {
  // output lineage of leaf;
  os << "{\"lineage\":[" << strjoin(lngs.begin(), lngs.end(), ",")
     << "]," << endl;

  // output the entropy
  os << "\"rank\":";
  aTaxa.rank->outRanksJson(os);
  os << "}" << endl;
};
