/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-11-27 09:59:06
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-09 21:51:04
 */

#include "getLineage.h"

int main(int argc, char *argv[]) {

  // get the input arguments
  Args myargs(argc, argv);

  // Initial the taxadb by the dump files
  TaxaDB taxdb(myargs.dbpath);

  // reset the rank abbrivation map
  TaxaRank *rank = TaxaRank::create();
  if (!myargs.rankfile.empty())
    rank->setRankByFile(myargs.rankfile);

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
    vector<string> hit;
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

    // output the result
    for (auto &it : hit) {
      ofs << it << endl;
    }
    if (myargs.outNonhit) {
      for (auto &it : nonhit)
        ofs << it << endl;
    }
  }
  ofs.close();
}

Args::Args(int argc, char **argv)
    : program(argv[0]), outfile("Lineage.list"), queryfile("name.list"),
      queryTaxID(0), queryName(""), outNonhit(true) {

  char ch;
  while ((ch = getopt(argc, argv, "i:d:o:r:I:N:Hqh")) != -1) {
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

  if (dbpath.empty()) {
    dbpath = "taxadb.gz";
    if (!fileExists(dbpath)) {
      dbpath = "./taxdump";
      if (!fileExists(dbpath)) {
        cerr << "\n** Error **: No database file/directory gived" << endl;
        usage();
      }
    }
  }
}

void Args::usage() {
  cerr << "\nProgram Usage: \n"
       << program << "\n"
       << " [ -i name.list ]    The query list file defalut: blank, \n"
       << "                     output all items of database \n"
       << " [ -d taxadb.gz ]    The dump of NCBI taxonomy database\n"
       << " [ -o Lineage.list ] Output file, default: stdout\n"
       << " [ -r <rankfile> ]   Rank mapping file\n"
       << " [ -q ]              Run command in quiet mode\n"
       << " [ -h ]              display this information\n"
       << endl;
  exit(1);
}

void twoColumn(TaxaDB &taxdb, const string &fname, vector<string> &hit,
               vector<string> &nonhit, int ncName, int ncTaxid) {
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
      hit.emplace_back(nmlist[i] + " " + res);
    }
  }
}

void oneColumn(TaxaDB &taxdb, const string &fname, vector<string> &hit,
               vector<string> &nonhit, int ncName) {
  vector<string> qlist;
  readlist(fname, qlist, ncName);
  // search the query in list by name
  for (auto &q : qlist) {
    string res = taxdb.search(q);
    if (res.empty()) {
      nonhit.push_back(q);
    } else {
      hit.emplace_back(q + " " + res);
    }
  }
};
