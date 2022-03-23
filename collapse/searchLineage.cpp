/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-23 09:25:45
 */

#include "searchLineage.h"

void searchLineage(int argc, char *argv[]) {

  // get the input arguments
  UpLngArgs myargs(argc, argv);

  // set rank
  TaxaRank *rank = TaxaRank::create();
  rank->initial(myargs.rankfile, myargs.outrank);

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
    : outfile("Lineage.csv"), format("") {

  program = argv[0];
  string infile("namelist.txt");

  char ch;
  while ((ch = getopt(argc, argv, "i:d:o:r:m:l:R:F:Jqh")) != -1) {
    switch (ch) {
    case 'i':
      infile = optarg;
      break;
    case 'd':
      taxadb = optarg;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'm':
      taxrev = optarg;
      break;
    case 'r':
      outrank = optarg;
      break;
    case 'R':
      rankfile = optarg;
      break;
    case 'l':
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

  // read the input name list
  smatch matchs;
  if (regex_search(infile, matchs, regex("([^:]+):([0-9]+)"))) {
    string fname = matchs[1].str();
    int ncol = stoi(matchs[2].str());
    readlist(fname, nmlist, ncol);
  } else {
    readlist(infile, nmlist, 1);
  }

  // set the output format
  if (format.empty()) {
    format = getsuffix(outfile);
  }

  // set the custom lineage file
  if (taxfile.empty()) {
    taxfile = "Lineage.txt";
    if (!fileExists(taxfile))
      taxfile = "Lineage.csv";
  }
  
  // for the default
  if (taxadb.empty()) {
    taxadb = "taxadb.gz";
    if (!fileExists(taxadb)) {
      taxadb = "taxdump.tar.gz";
      if (!fileExists(taxadb)) {
        taxadb = "taxdump";
      }
    }
  }
}

void UpLngArgs::usage() {
  cerr << "\nProgram Usage: \n\n"
       << program << "\n"
       << " [ -i namelist.txt ]    Input name list, ':N' after the file name\n"
       << "                        select the N column of the file\n"
       << "                        default: first column of namelist.txt\n"
       << " [ -o Lineage.csv ]     Output lineage file, default: lineage.csv\n"
       << " [ -m <Revision.txt> ]  Lineage revise file for batch edit,\n"
       << "                        default: None\n"
       << " [ -l Lineage.txt ]     Lineage file for leafs of tree, \n"
       << "                        default: Lineage.txt\n"
       << " [ -d taxadb.gz ]       Taxa database file or directory,\n"
       << "                        default: taxadb.gz or taxdump.tar.gz\n"
       << " [ -R <None> ]          List file for rank names and abbreviations,\n"
       << "                        default: use the setting of program\n"
       << " [ -r <DKPCOFGS> ]      Set output taxon rank by abbreviations,\n"
       << "                        default: according to source\n"
       << " [ -q ]                 Run command in quiet mode\n"
       << " [ -h ]                 Display this information\n"
       << endl;
  exit(1);
}