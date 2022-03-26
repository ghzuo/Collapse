/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-26 13:59:26
 */

#include "getLeafName.h"

void getLeafName(int argc, char *argv[]) {

  // get the name of file
  string infile("Tree.nwk");
  string outfile("namelist.txt");
  string program(argv[0]);

  char ch;
  while ((ch = getopt(argc, argv, "i:o:h")) != -1) {
    switch (ch) {
    case 'i':
      infile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'h':
      leafUsage(program);
    case '?':
      leafUsage(program);
    }
  }

  Node *aTree = new Node;
  aTree->innwk(infile);
  vector<Node *> allLeafs;
  aTree->getLeafs(allLeafs);

  ofstream is(outfile);
  for (auto &nd : allLeafs) {
    is << nd->name << endl;
  }
  is.close();
}

void leafUsage(string &program) {
  cerr << "\nProgram Usage: \n"
       << program << "\n"
       << " [ -i Tree.nwk ]      Input tree file, default: Tree.nwk\n"
       << " [ -o namelist.txt ]  Output name list file, \n"
       << "                      default: namelist.txt\n"
       << " [ -q ]               Run command in quiet mode\n"
       << " [ -h ]               Display this information\n"
       << endl;
  exit(1);
}