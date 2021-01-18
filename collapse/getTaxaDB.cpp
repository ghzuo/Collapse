/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-12-08 20:01:45
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2021-01-18 14:49:45
 */

#include "getTaxaDB.h"

void mkTaxaDB(int argc, char *argv[]) {

  // get the name of file
  string dumps("taxdump.tar.gz");
  string outfile("taxadb.gz");
  string program(argv[0]);

  char ch;
  while ((ch = getopt(argc, argv, "i:o:h")) != -1) {
    switch (ch) {
    case 'i':
      dumps = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'h':
      mkdbUsage(program);
    case '?':
      mkdbUsage(program);
    }
  }

  // Initial the taxadb by the dump files
  TaxaDB taxdb(dumps);

  // output the gz database if queryfile empty
  taxdb.writeTable(outfile);
}

void mkdbUsage(string &program) {
  cerr
      << "\nProgram Usage: \n"
      << program << "\n"
      << " [ -d taxdump.tar.gz ]  NCBI taxon dumpfile directory, default: taxdump\n"
      << " [ -o taxadb.gz ]       Packaged taxon database, default: taxadb.gz\n"
      << " [ -q ]                 Run command in quiet mode\n"
      << " [ -h ]                 Display this information\n"
      << endl;
  exit(1);
}