/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-16 12:26:16
 */

#include "updateLineage.h"

void updateLineage(int argc, char *argv[]) {

  // get the input arguments
  UpLngArgs myargs(argc, argv);

  // set rank
  TaxaRank *rank = TaxaRank::create();
  rank->initial(myargs.abfile, myargs.abtype);

  // set the lineage
  LngData lngs(myargs.taxadb, myargs.taxfile, myargs.taxrev);
  lngs.getLineage(myargs.nmlist);

  // output data
  lngs.output(myargs.outfile, myargs.format);
}
/****************************************************************************
 ******************************  End main program ***************************
 ****************************************************************************/

UpLngArgs::UpLngArgs(int argc, char **argv)
    : taxrev("Lineage.rev"), taxfile("Lineage.list"), taxadb("taxadb.gz"),
      outfile("Lineage.csv"), format("") {

  program = argv[0];
  string infile("name.list");

  char ch;
  while ((ch = getopt(argc, argv, "i:D:o:r:t:T:R:F:Jqh")) != -1) {
    switch (ch) {
    case 'i':
      infile = optarg;
      break;
    case 'D':
      taxadb = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'r':
      taxrev = optarg;
      break;
    case 't':
      abtype = optarg;
      break;
    case 'R':
      abfile = optarg;
      break;
    case 'T':
      taxfile = optarg;
      break;
    case 'F':
      format = optarg;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'J':
      theJson.format = true;
      break;
    case 'h':
      usage();
    case '?':
      usage();
    }
  }

  smatch matchs;
  if (regex_search(infile, matchs, regex("([^:]+):([0-9]+)"))) {
    string fname = matchs[1].str();
    int ncol = stoi(matchs[2].str());
    readlist(fname, nmlist, ncol);
  } else {
    readlist(infile, nmlist, 1);
  }

  if (format.empty()) {
    format = getsuffix(outfile);
  }
}

void UpLngArgs::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -i name.list ]     Input name list, ':N' after the file name\n"
       << "                      select the N column of the file\n"
       << "                      default: first column of name.list\n"
       << " [ -o Lineage.csv ]   Output prefix name: default: collapsed\n"
       << " [ -r Lineage.rev ]   Lineage revise file for batch edit,\n"
       << "                      default: Lineage.rev\n"
       << " [ -T Lineage.list ]  Lineage file for leafs of tree, \n"
       << "                      default: Lineage.list\n"
       << " [ -D taxadb.gz ]     Taxa database file or directory,\n"
       << "                      default: taxadb.gz or taxdump/\n"
       << " [ -R <None> ]        Abbravition list for taxon rank name,\n"
       << "                      default: by program\n"
       << " [ -t DKPCOFGS ]      Abbreviation of output taxon rank,\n"
       << "                      default: DKPCOFGS\n"
       << " [ -q ]               Run command in quiet mode\n"
       << " [ -h ]               Display this information\n"
       << endl;
  exit(1);
}