/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-11-27 09:59:06
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-06 12:38:29
 */

#ifndef GETLINEAGE_H
#define GETLINEAGE_H

#include "taxadb.h"
#include <string>

using namespace std;

struct Args {
  string program;
  string indir;
  string namefile;
  string nodefile;
  string outfile;
  string rankfile;
  vector<string> qlist;

  Args(int, char **);
  void usage(string &);
};

regex regTaxID{"taxid([0-9]+)"};
regex regBinomial{"^([A-Z][a-z]+_[a-z]+)[_.]"};
regex regUnitName{"^([A-Z][a-z]+)"};

string searchLineage(const string &);

#endif // !GETLINEAGE_H