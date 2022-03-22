/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-22 19:22:34
 */

#ifndef GETLINEAGE_H
#define GETLINEAGE_H

#include "taxadb.h"
#include <regex>
#include <string>

using namespace std;

struct LngArgs {
  string program;
  string dbpath;
  string outfile;
  string rankfile;
  string outrank;
  string queryfile;
  size_t queryTaxID;
  string queryName;
  bool outNonhit;

  LngArgs(int, char **);
  void usage();
};

void queryLineage(int, char **);

void twoColumn(TaxaDB &, const string &, vector<pair<string,string>> &, vector<string> &,
               int ncName = 1, int ncTaxid = 2);
void oneColumn(TaxaDB &, const string &, vector<pair<string,string>> &, vector<string> &,
               int ncName = 1);

#endif // !GETLINEAGE_H