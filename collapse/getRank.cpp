/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-22 15:10:34
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-26 23:38:00
 */

#include "getRank.h"
void getRank(int argc, char *argv[]) {

  // get the name of file
  string outfile("ranklist.txt");
  string program(argv[0]);
  bool allrank(false);

  char ch;
  while ((ch = getopt(argc, argv, "o:aqh")) != -1) {
    switch (ch) {
    case 'o':
      outfile = optarg;
      break;
    case 'a':
      allrank = true;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'h':
      rankUsage(program);
    case '?':
      rankUsage(program);
    }
  }

  // get the rank abbreviation map
  TaxaRank *rank = TaxaRank::create();
  ofstream ofs(outfile);

  if (allrank) {
    for (auto &rk : rank->rankmap) {
      ofs << rk.first << "\t" << rk.second << endl;
    }
  } else {
    vector<string> ranklist = {"kingdom", "phylum", "class",  "order",
                               "family",  "genus",  "species"};
    vector<string> prefix{"super", "", "sub"};

    for (auto nm : ranklist) {
      for (auto pf : prefix) {
        string name = pf + nm;
        auto iter = rank->rankmap.find(name);
        if (iter != rank->rankmap.end()) {
          ofs << iter->first << "\t" << iter->second << endl;
        }
      }
    }
  }
  ofs.close();
}

void rankUsage(string &program) {
  cerr << "\nProgram Usage: \n"
       << program << "\n"
       << " [ -o ranklist.txt ]  Output name list file,\n"
       << "                      default: ranklist.txt\n"
       << " [ -a ]               Output all inbuilt ranks\n"
       << " [ -q ]               Run command in quiet mode\n"
       << " [ -h ]               Display this information\n"
       << endl;
  exit(1);
}