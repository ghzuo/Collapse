/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-06-01 16:23:51
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-06-01 16:39:13
 */

#include "rooting.h"

void rooting(int argc, char *argv[]) {

  // get the input arguments
  RootingArgs myargs(argc, argv);

  /************************************************************************
   ******* read the tree and rooting by topology and branch ***************/
  Node *aTree = new Node;
  aTree->innwk(myargs.infile);

  // rooting the tree by outgroup and branch length
  if (aTree->children.size() == 2) {
    theInfo("This tree is a rooted tree, keep as it is");
  } else if (!myargs.outgrp.empty()) {
    // root the unrooted tree by input outgroup name
    Node *result = aTree->rootingByOutgrp(myargs.outgrp);
    if (result == NULL) {
      theInfo("Rooting the tree by branch length");
      aTree = aTree->rootingByLength(myargs.rootMeth);
    } else {
      aTree = result;
    }
  } else {
    aTree = aTree->rootingByLength(myargs.rootMeth);
  }

  // output tree
  aTree->outnwk(myargs.outfile);
}

RootingArgs::RootingArgs(int argc, char **argv)
    : infile("Tree.nwk"), outfile("Rooted.nwk"),outgrp(""), rootMeth("mv"){

  program = argv[0];

  char ch;
  while ((ch = getopt(argc, argv,
                      "i::o:O::m:u:qh")) != -1) {
    switch (ch) {
    case 'i':
      infile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'O':
      outgrp = optarg;
      break;
    case 'm':
      rootMeth = toLower(optarg);
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
}

void RootingArgs::usage() {
  cerr
      << "\nProgram Usage: \n\n"
      << program << "\n"
      << " [ -i Tree.nwk ]     Input newick tree, default: Tree.nwk\n"
      << " [ -o Rooted.nwk ]   Output rooted newick tree, default: Rooted.nwk \n"
      << " [ -m mv ]           Set rooting method: mv, mad, mp, pmr, or md\n"
      << "                     default: mv\n"
      << " [ -O <Outgroup> ]   Rooting phylogenetic tree by outgroup.\n"
      << "                     default: None, by mv method\n"
      << " [ -P ]              Output prediction for undefined leafs\n"
      << " [ -q ]              Run command in quiet mode\n"
      << " [ -h ]              Display this information\n"
      << endl;
  exit(1);
}
