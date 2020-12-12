/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-09-01 12:54:33
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-12 14:44:36
 */

#include "collapse.h"

int main(int argc, char *argv[]) {

  // get the input arguments
  Args myargs(argc, argv);

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
  output(aTaxa, aTree, lngs, myargs);
}
/****************************************************************************
 ******************************  End main program ***************************
 ****************************************************************************/

Args::Args(int argc, char **argv)
    : infile(""), taxrev(""), outgrp(""), taxfile(""), lngfile(""),
      forWeb(false), predict(true) {

  program = argv[0];
  string outname("collapsed");
  string supdir("./");
  string kstr;

  char ch;
  while ((ch = getopt(argc, argv, "i:d:D:o:r:t:s:T:L:R:O:WPqh")) != -1) {
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
    case 'L':
      lngfile = optarg;
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

  // for the output newick file name
  treeSuff = ".tree";
  if (infile.find(treeSuff) + treeSuff.size() == infile.size())
    treeSuff = ".nwk";

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

void Args::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -d ./ ]            The work directory, default: ./\n"
       << " [ -i Tree.nwk ]      Input newick tree file, default: Tree.nwk\n"
       << " [ -o collpsed ]      Output prefix name: default: collapsed\n"
       << " [ -r Lineage.rev ]   Lineage revise file for batch edit,\n"
       << "                      default: Lineage.rev\n"
       << " [ -T Lineage.list ]  Lineage file for leafs of tree, \n"
       << "                      default: Lineage.list\n"
       << " [ -D taxadb.gz ]     Taxa database file or directory,\n"
       << "                      default: taxadb.gz or taxdump/\n"
       << " [ -L <None> ]        Output lineage of all leafs:\n"
       << "                      default: don't output lineages\n"
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
void output(Taxa &aTaxa, Node *aTree, vector<Lineage> &lngs, Args &myargs) {

  /// output the monophyly
  aTaxa.outStatitics(myargs.outPref + ".unit");

  /// output for webserver
  if (myargs.forWeb) {
    //... check the upload genomes
    aTree->checkUploaded();

    /// output the json file
    aTree->outjson(myargs.outPref + ".json");

    /// output the taxonomy list with relationship
    aTaxa.outTax(myargs.outPref + ".list");
  } else {
    /// output the newick file
    aTree->outnwk(myargs.outPref + myargs.treeSuff);

    /// output the entropy
    aTaxa.outEntropy(myargs.outPref + ".entropy");
  }

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

  // output lineage of leaf;
  if (!myargs.lngfile.empty()) {
    ofstream leaf(myargs.lngfile);
    for (auto &lng : lngs) {
      leaf << lng.name << endl;
    }
    leaf.close();
  }
};
