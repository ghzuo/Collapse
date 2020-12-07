/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-11-27 09:59:06
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-06 15:05:02
 */

#include "getLineage.h"

int main(int argc, char *argv[]) {

  // get the input arguments
  Args myargs(argc, argv);

  // Initial the taxadb by the dump files
  TaxaDB tdb(myargs.nodefile, myargs.namefile);

  // reset the rank abbrivation map
  if (!myargs.rankfile.empty())
    tdb.resetRankMap(myargs.rankfile);


  ofstream ofs(myargs.outfile);
  if (myargs.qlist.empty()) {
    tdb.exportLineage(ofs);
  } else {
    for (auto &q : myargs.qlist) {
      ofs << q << "\t" << tdb.search(q) << endl;
    }
  }
  ofs.close();
}

Args::Args(int argc, char **argv)
    : indir("./"), outfile("Lineage.list"), program(argv[0]) {

  string queryfile;

  char ch;
  while ((ch = getopt(argc, argv, "i:n:d:o:r:qh")) != -1) {
    switch (ch) {
    case 'd':
      indir = optarg;
      addsuffix(indir, '/');
      break;
    case 'i':
      queryfile = optarg;
      break;
    case 'n':
      nodefile = optarg;
      break;
    case 'm':
      namefile = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'r':
      rankfile = optarg;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'h':
      usage(program);
    case '?':
      usage(program);
    }
  }

  if (nodefile.empty())
    nodefile = indir + "nodes.dmp";
  if (namefile.empty())
    namefile = indir + "names.dmp";

  if (!queryfile.empty())
    readlist(queryfile, qlist);
}

void Args::usage(string &program) {
  cerr << "\nProgram Usage: \n"
       << program << "\n"
       << " [ -d ./ ]           the directory of the dump files directory\n"
       << " [ -i qfile ]        the query list file"
       << " [ -n nodes.dmp ]    the file of nodes.dmp\n"
       << " [ -m names.dmp ]    the file of names.dmp\n"
       << " [ -o <outfile>]     output file, default: stdout\n"
       << " [ -r rankfile ]     Rank mapping file\n"
       << " [ -q ]              Run command in quiet mode\n"
       << " [ -h ]              display this information\n"
       << endl;
  exit(1);
}
