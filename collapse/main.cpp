/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-26 14:08:23
 */

#include "collapse.h"
#include "getLeafName.h"
#include "getRank.h"
#include "getTaxaDB.h"
#include "queryLineage.h"
#include "searchLineage.h"
using namespace std;

void usage(string &program) {
  cerr << "\nProgram Usage: \n"
       << program << " Token [options] \n"
       << " Available Tokens:\n"
       << "   run       All-in-one command: search lineage of leavies,\n"
       << "             annotate phylogenetic tree, and do comparation.\n"
       << "   leaf      Obtain species name list of phylogenetic tree\n"
       << "   cache     Make NCBI database cache from taxdump.tar.gz\n"
       << "   rank      Output taxon rank names and abbreviations\n"
       << "   query     Query lineage from local NCBI taxonomy database\n"
       << "   search    Search lineage from lineage files and NCBI taxonomy\n"
       << "             database, and revised by the revision file\n"
       << "   help      Provide the help information for <Token>\n"
       << " [ -h ]      Display this information\n"
       << endl;
  exit(1);
}

int main(int argc, char *argv[]) {

  // get the program and task
  string program(argv[0]);
  string task;
  if (argc > 1)
    task = argv[1];

  if (task.empty() || task[0] == '-') {
    if (task.compare("-h") == 0) {
      usage(program);
    } else {
      collapse(argc, argv);
    }
  } else {
    // set the argc and argv for task
    *(argv[1] - 1) = ' ';
    argv[1] = argv[0];
    argv = &argv[1];
    argc -= 1;

    // for help task
    if (task.compare("help") == 0) {
      if (argc > 1) {
        task = argv[1];
        *(argv[1] - 1) = ' ';
        argv[1] += (strlen(argv[1]) + 1);
        strcpy(argv[1], "-h");
        argc = 2;
      } else {
        usage(program);
      }
    }

    // run the task
    if (task.compare("run") == 0) {
      collapse(argc, argv);
    } else if (task.compare("query") == 0) {
      queryLineage(argc, argv);
    } else if (task.compare("cache") == 0) {
      mkTaxaDB(argc, argv);
    } else if (task.compare("leaf") == 0) {
      getLeafName(argc, argv);
    } else if (task.compare("search") == 0) {
      searchLineage(argc, argv);
    } else if (task.compare("rank") == 0) {
      getRank(argc, argv);
    } else if (task.compare("help") == 0) {
      usage(program);
    } else {
      cerr << "\nUnknown Task: " << task << endl;
      usage(program);
    }
  }
}
