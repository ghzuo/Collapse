/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-11-27 09:59:06
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2021-01-18 14:30:55
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
  string queryfile;
  size_t queryTaxID;
  string queryName;
  bool outNonhit;

  LngArgs(int, char **);
  void usage();
};

void getlineage(int, char **);

void twoColumn(TaxaDB &, const string &, vector<string> &, vector<string> &,
               int ncName = 1, int ncTaxid = 2);
void oneColumn(TaxaDB &, const string &, vector<string> &, vector<string> &,
               int ncName = 1);

#endif // !GETLINEAGE_H