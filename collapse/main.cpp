/*
 * Copyright (c) 2021  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2021-01-18 10:31:12
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2021-01-18 16:01:07
 */

#include "collapse.h"
#include "getLeafName.h"
#include "getLineage.h"
#include "getTaxaDB.h"
using namespace std;

void usage(string &program) {
  cerr << "\nProgram Usage: \n"
       << program << " Task [options] \n"
       << " Available Task:\n"
       << "   run       Compare phylogenetic tree with taxonomy system\n"
       << "   index     Index local NCBI taxonomy by dump file\n"
       << "   query     Query lineage from the taxonomy database\n"
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
  } else {
    cerr << "\nUnknown Task: " << task << endl;
    usage(program);
  }
}