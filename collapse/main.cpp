/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-23 11:54:11
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
       << program << " Task [options] \n"
       << " Available Task:\n"
       << "   run       Annotate phylogenetic tree with taxonomy system\n"
       << "   cache     Make NCBI database cache from taxdump.tar.gz\n"
       << "   query     Query lineage from local NCBI taxonomy database\n"
       << "   search    Search lineage from lineage files and NCBI taxonomy \n"
       << "             database, and revised by the revision file\n"
       << "   leaf      Obtain the species name list of phylogenetic tree\n"
       << "   rank      Output the default rank names and abbreviations\n"
       << "   help      Privode the help information for <Task>\n"
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
    char **subArgv = &argv[1];
    int subArgc = argc - 1;

    // add help task
    if (task.compare("help") == 0) {
      if (subArgc > 1) {
        task = subArgv[1];
        subArgc = 2;
        *(argv[2] - 1) = ' ';
        int len = strlen(argv[2]);
        strcat(argv[2], " -h");
        argv[2][len] = '\0';
        argv[3] = &(argv[2][len + 1]);
        argv[2] = argv[0];
        subArgv = &(argv[2]);
      } else {
        usage(program);
      }
    }

    // run the task
    if (task.compare("run") == 0) {
      collapse(subArgc, subArgv);
    } else if (task.compare("query") == 0) {
      queryLineage(subArgc, subArgv);
    } else if (task.compare("cache") == 0) {
      mkTaxaDB(subArgc, subArgv);
    } else if (task.compare("leaf") == 0) {
      getLeafName(subArgc, subArgv);
    } else if (task.compare("search") == 0) {
      searchLineage(subArgc, subArgv);
    } else if (task.compare("rank") == 0) {
      getRank(subArgc, subArgv);
    } else if (task.compare("help") == 0) {
      usage(program);
    } else {
      cerr << "\nUnknown Task: " << task << endl;
      usage(program);
    }
  }
}