/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-09-01 12:54:33
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-07 16:52:33
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
  LineageHandle lngHandle(myargs.taxfile, myargs.taxrev, myargs.abfile,
                          myargs.abtype);

  // get the leafs lineage
  vector<Node *> allLeafs;
  aTree->getLeafs(allLeafs);
  vector<Lineage> lngs;
  for (auto nd : allLeafs) {
    lngs.emplace_back(nd->name);
  }

  // set lineage of leafs
  lngHandle.getLineage(lngs);

  /*************************************************************************
   *********** set lineage of Nodes and rooting unrooted tree **************/
  // set lineage of leaf
  for (int i = 0; i < allLeafs.size(); ++i) {
    allLeafs[i]->setOneLeaf(lngs[i]);
  }

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

  /**************************************************************************
   ************ get the statistic of Taxonomy and annotate tree *************/
  TaxSys::setDivision(lngHandle.outrank);

  // get the statics of the two taxa system (defined and undefined)
  Taxa aTaxa(lngs);

  // annotate taxonomy of nodes
  aTaxa.annotate(aTree);

  /***************************************************************************
   *********  output data ****************************************************/
  output(aTaxa, aTree, myargs);
}
/****************************************************************************
 ******************************  End main program ***************************
 ****************************************************************************/

Args::Args(int argc, char **argv)
    : infile(""), taxrev(""), outgrp(""), taxfile(""), forWeb(false),
      predict(true) {

  program = argv[0];
  string outname("collapsed");
  string supdir("./");
  string kstr;

  char ch;
  while ((ch = getopt(argc, argv, "i:d:o:r:t:s:T:L:O:WPqh")) != -1) {
    switch (ch) {
    case 'i':
      infile = optarg;
      break;
    case 'd':
      supdir = optarg;
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
    case 'L':
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
}

void Args::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -d ./ ]           The work directory, default: ./\n"
       << " [ -i Tree.nwk ]     Input cv file list, default: Tree.nwk\n"
       << " [ -o collpsed ]     Output dist: default: collapsed\n"
       << " [ -r Lineage.rev ]  Lineage revise file: default: Lineage.rev\n"
       << " [ -L abfile ]       Unusual abbravition list for taxon level name\n"
       << " [ -t DKPCOFGS ]     Abbreviation of taxonomy level, default: "
          "DKPCOFGS\n"
       << " [ -O <Outgroup> ]   Set the outgroup for the unroot tree.\n"
       << " [ -T taxlist ]      Taxonomy map files, default: None\n"
       << " [ -q ]              Run command in quiet mode\n"
       << " [ -h ]              Display this information\n"
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
void output(Taxa &aTaxa, Node *aTree, Args &myargs) {

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
};
