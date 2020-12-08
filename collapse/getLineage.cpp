/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-11-27 09:59:06
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-08 10:24:26
 */

#include "getLineage.h"

int main(int argc, char *argv[]) {

  // get the input arguments
  Args myargs(argc, argv);

  // Initial the taxadb by the dump files
  TaxaDB taxdb(myargs.dbpath);

  // reset the rank abbrivation map
  if (!myargs.rankfile.empty())
    taxdb.resetRankMap(myargs.rankfile);

  // do search lineage
  if (myargs.queryfile.empty()) {
    taxdb.writeTable(myargs.outfile);
  } else {
    ofstream ofs(myargs.outfile);
    smatch matchs;
    regex_search(myargs.queryfile, matchs, regex("([^:]+):([0-9]+),([0-9]+)"));
    if (!matchs.empty()) {
      string fname = matchs[1].str();
      int ncName = stoi(matchs[2].str());
      int ncTaxid = stoi(matchs[3].str());
      multiColumns(taxdb, fname, ofs, ncName, ncTaxid);
    } else {
      regex_search(myargs.queryfile, matchs, regex("([^:]+):([0-9]+)"));
      if (!matchs.empty()) {
        string fname = matchs[1].str();
        int ncName = stoi(matchs[2].str());
        oneColumn(taxdb, fname, ofs, ncName);
      } else if (nColumns(myargs.queryfile) > 1) {
        multiColumns(taxdb, myargs.queryfile, ofs);
      } else {
        oneColumn(taxdb, myargs.queryfile, ofs);
      }
    }
    ofs.close();
  }
}

Args::Args(int argc, char **argv)
    : dbpath("./taxdump"), outfile("Lineage.list"), program(argv[0]) {

  char ch;
  while ((ch = getopt(argc, argv, "i:d:o:r:qh")) != -1) {
    switch (ch) {
    case 'd':
      dbpath = optarg;
      break;
    case 'i':
      queryfile = optarg;
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
}

void Args::usage(string &program) {
  cerr << "\nProgram Usage: \n"
       << program << "\n"
       << " [ -i qfile ]        The query list file defalut: blank, \n"
       << "                     output all items of database \n"
       << " [ -d ./taxdump ]    The dump of NCBI taxonomy database\n"
       << " [ -o Lineage.list ] Output file, default: stdout\n"
       << " [ -r <rankfile> ]   Rank mapping file\n"
       << " [ -q ]              Run command in quiet mode\n"
       << " [ -h ]              display this information\n"
       << endl;
  exit(1);
}

void multiColumns(TaxaDB &taxdb, const string &fname, ostream &os, int ncName,
                  int ncTaxid) {
  // for two column file
  vector<string> nmlist;
  vector<size_t> tidlist;
  readlist(fname, nmlist, ncName);
  readlist(fname, tidlist, ncTaxid);
  for (size_t i = 0; i < nmlist.size(); ++i) {
    os << nmlist[i] << "\t" << taxdb.getLineage(tidlist[i]).lineage << endl;
  }
}

void oneColumn(TaxaDB &taxdb, const string &fname, ostream &os, int ncName) {
  vector<string> qlist;
  readlist(fname, qlist, ncName);
  for (auto &q : qlist) {
    os << q << "\t" << taxdb.search(q) << endl;
  }
};
