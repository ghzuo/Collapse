/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-26 13:44:15
 */

#include "getTaxaDB.h"

void mkTaxaDB(int argc, char *argv[]) {

  // get the name of file
  string dumps("taxdump.tar.gz");
  string outfile("taxadb.gz");
  string outdir("");
  string program(argv[0]);

  char ch;
  while ((ch = getopt(argc, argv, "i:o:J:h")) != -1) {
    switch (ch) {
    case 'i':
      dumps = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'J':
      outdir = optarg;
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
  if (outdir.empty()) {
    taxdb.writeTable(outfile);
  } else {
    theJson.format = true;
    taxdb.outJsons(outdir);
  }
}

void mkdbUsage(string &program) {
  cerr
      << "\nProgram Usage: \n"
      << program << "\n"
      << " [ -d taxdump.tar.gz ]  NCBI taxon dump, default: taxdump.tar.gz\n"
      << " [ -o taxadb.gz ]       Packaged taxon database, default: taxadb.gz\n"
      << " [ -q ]                 Run command in quiet mode\n"
      << " [ -h ]                 Display this information\n"
      << endl;
  exit(1);
}