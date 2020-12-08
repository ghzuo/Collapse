/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-11-27 09:59:06
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-07 20:24:06
 */

#ifndef GETLINEAGE_H
#define GETLINEAGE_H

#include "taxadb.h"
#include <string>
#include <regex>

using namespace std;

struct Args {
  string program;
  string dbpath;
  string outfile;
  string rankfile;
  string queryfile;

  Args(int, char **);
  void usage(string &);
};

void multiColumns(TaxaDB &, const string &, ostream &, int ncName = 1,
                  int ncTaxid = 2);
void oneColumn(TaxaDB &, const string &, ostream &, int ncName = 1);

#endif // !GETLINEAGE_H