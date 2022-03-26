/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-26 14:06:41
 */

#include "queryLineage.h"

void queryLineage(int argc, char *argv[]) {

  // get the input arguments
  LngArgs myargs(argc, argv);

  // Initial the taxadb by the dump files
  TaxaDB taxdb(myargs.dbpath);

  // reset the rank abbreviation map
  TaxaRank *rank = TaxaRank::create();
  if (!myargs.rankfile.empty())
    rank->setRankByFile(myargs.rankfile);
  if (!myargs.outrank.empty())
    rank->setOutRank(myargs.outrank);

  // do the search
  ofstream ofs(myargs.outfile);
  if (myargs.queryTaxID != 0 || !myargs.queryName.empty()) {
    // search one item
    vector<TaxonNode> tnlist;
    if (myargs.queryTaxID != 0) {
      taxdb.query(myargs.queryTaxID, tnlist);
    } else {
      taxdb.query(myargs.queryName, tnlist);
    }

    // output result
    for (auto &tn : tnlist) {
      ofs << tn << endl;
    }

  } else {
    // search the queries in query file
    smatch matchs;
    vector<pair<string, string>> hit;
    vector<string> nonhit;
    if (regex_search(myargs.queryfile, matchs,
                     regex("([^:]+):([0-9]+),([0-9]+)"))) {
      string fname = matchs[1].str();
      int ncName = stoi(matchs[2].str());
      int ncTaxid = stoi(matchs[3].str());
      twoColumn(taxdb, fname, hit, nonhit, ncName, ncTaxid);
    } else if (regex_search(myargs.queryfile, matchs,
                            regex("([^:]+):([0-9]+)"))) {
      string fname = matchs[1].str();
      int ncName = stoi(matchs[2].str());
      oneColumn(taxdb, fname, hit, nonhit, ncName);
    } else if (nColumns(myargs.queryfile) > 1) {
      twoColumn(taxdb, myargs.queryfile, hit, nonhit);
    } else {
      oneColumn(taxdb, myargs.queryfile, hit, nonhit);
    }

    // format the output lineage string
    if (!myargs.outrank.empty()) {
      for (auto &it : hit) {
        rank->format(it.second);
      }
    }

    // output the result
    for (auto &it : hit) {
      ofs << it.first << " " << it.second << endl;
    }

    // output the nonhit
    if (myargs.outNonhit) {
      for (auto &it : nonhit)
        ofs << it << endl;
    }
  }
  ofs.close();
}

LngArgs::LngArgs(int argc, char **argv)
    : program(argv[0]), outfile("Lineage.txt"), queryfile("namelist.txt"),
      queryTaxID(0), queryName(""), outNonhit(true) {

  char ch;
  while ((ch = getopt(argc, argv, "i:d:o:R:r:I:N:Hqh")) != -1) {
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
    case 'R':
      rankfile = optarg;
      break;
    case 'r':
      outrank = optarg;
      break;
    case 'I':
      queryTaxID = stoi(optarg);
      break;
    case 'N':
      queryName = optarg;
      break;
    case 'H':
      outNonhit = false;
      break;
    case 'q':
      theInfo.quiet = true;
      break;
    case 'h':
      usage();
    case '?':
      usage();
    }
  }

  // check the input file
  if (!fileExists(queryfile)) {
    cerr << "\n** Error **: Cannot find the input file : " << queryfile << endl;
    usage();
  }

  // set the default NCBI Taxonomy database file/folder
  if (dbpath.empty()) {
    dbpath = "taxadb.gz";
    if (!fileExists(dbpath)) {
      dbpath = "taxdump.tar.gz";
      if (!fileExists(dbpath)) {
        dbpath = "./taxdump";
        if (!fileExists(dbpath)) {
          cerr << "\n** Error **: No database file/directory gived" << endl;
          usage();
        }
      }
    }
  }
}

void LngArgs::usage() {
  cerr << "\nProgram Usage: \n"
       << program << "\n"
       << " [ -I <Taxon ID> ]     Query a taxon id\n"
       << " [ -N <Taxon Name> ]   Query a taxon name\n"
       << " [ -i namelist.txt ]   The query list file,\n"
       << "                       defalut: namelist.txt\n"
       << " [ -d taxadb.gz ]      The NCBI taxonomy data file\n"
       << "                       default: taxadb.gz or taxdump.tar.gz\n"
       << " [ -o Lineage.txt ]    Output lineage file, \n"
       << "                       default: Lineage.txt\n"
       << " [ -R <None> ]         List file for rank names and abbreviations,\n"
       << "                       default: set by program\n"
       << " [ -r DKPCOFGS ]       Set output taxon rank by abbreviations,\n"
       << "                       default: set by program\n"
       << " [ -H ]                Don't output missing items\n"
       << " [ -q ]                Run command in quiet mode\n"
       << " [ -h ]                Display this information\n"
       << endl;
  exit(1);
}

void twoColumn(TaxaDB &taxdb, const string &fname,
               vector<pair<string, string>> &hit, vector<string> &nonhit,
               int ncName, int ncTaxid) {
  // for two column file
  vector<string> nmlist;
  vector<size_t> tidlist;
  readlist(fname, nmlist, ncName);
  readlist(fname, tidlist, ncTaxid);
  // search the query in list by taxid column
  for (size_t i = 0; i < nmlist.size(); ++i) {
    string res = taxdb.search(tidlist[i]);
    if (res.empty()) {
      nonhit.push_back(nmlist[i]);
    } else {
      hit.emplace_back(nmlist[i], res);
    }
  }
}

void oneColumn(TaxaDB &taxdb, const string &fname,
               vector<pair<string, string>> &hit, vector<string> &nonhit,
               int ncName) {
  vector<string> qlist;
  readlist(fname, qlist, ncName);
  // search the query in list by name
  for (auto &q : qlist) {
    string res = taxdb.search(q);
    if (res.empty()) {
      nonhit.push_back(q);
    } else {
      hit.emplace_back(q, res);
    }
  }
};
