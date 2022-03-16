/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-16 12:29:13
 */

#include "collapse.h"
#include "getLeafName.h"
#include "getLineage.h"
#include "getTaxaDB.h"
#include "updateLineage.h"
using namespace std;

void usage(string &program) {
  cerr << "\nProgram Usage: \n"
       << program << " Task [options] \n"
       << " Available Task:\n"
       << "   run       Compare phylogenetic tree with taxonomy system\n"
       << "   index     Index local NCBI taxonomy by dump file\n"
       << "   query     Query lineage from the taxonomy database\n"
       << "   search    Search lineage from taxa file and taxonomy database\n"
       << "   leaf      Obtain the name list of phylogenetic tree\n"
       << " [ -h ]      Display this information\n"
       << endl;
  exit(1);
}

int main(int argc, char *argv[]) {

  // get the program and task
  string program(argv[0]);
  string task(argv[1]);

  // set the argc and argv for task
  *(argv[1] -1) = ' ';
  argv[1] = argv[0];
  char **subArgv = &argv[1];
  int subArgc = argc - 1;

  // run the task
  if (task.compare("run") == 0) {
    collapse(subArgc, subArgv);
  } else if (task.compare("query") == 0) {
    getlineage(subArgc, subArgv);
  } else if (task.compare("index") == 0) {
    mkTaxaDB(subArgc, subArgv);
  } else if (task.compare("leaf") == 0) {
    getLeafName(subArgc, subArgv);
  } else if (task.compare("search") == 0) {
    updateLineage(subArgc, subArgv);
  } else {
    cerr << "\nUnknown Task: " << task << endl;
    usage(program);
  }
}