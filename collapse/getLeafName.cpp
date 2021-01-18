/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-12-08 20:01:45
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2021-01-18 14:57:15
 */

#include "getLeafName.h"

void getLeafName(int argc, char *argv[]) {

  // get the name of file
  string infile("Tree.nwk");
  string outfile("name.list");
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
       << " [ -i Tree.nwk ]   Input tree file, default: Tree.nwk\n"
       << " [ -o name.list ]  Output name list, default: name.list\n"
       << " [ -q ]            Run command in quiet mode\n"
       << " [ -h ]            Display this information\n"
       << endl;
  exit(1);
}