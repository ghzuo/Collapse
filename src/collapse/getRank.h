/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-22 14:58:38
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-22 16:33:59
 */

#ifndef GETRANK_H
#define GETRANK_H

#include "taxadb.h"
#include <regex>
#include <string>

using namespace std;

void rankUsage(string &program);
void getRank(int, char**);

#endif //!.GETRANK_H