/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-12-12 10:07:46
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-12 12:16:20
 */

#ifndef FILEOPT_H
#define FILEOPT_H

#include <iostream>
#include <string>
#include <zlib.h>
using namespace std;

/********************************************************************************
 * @brief options for tar and zlib
 *
 ********************************************************************************/
#define  RECORDSIZE  512
struct TarHeader {
  char name[100];
  char mode[8];
  char uid[8];
  char gid[8];
  char size[12];
  char mtime[12];
  char chksum[8];
  char linkflag;
  char linkname[100];
  char magic[8];
  char uname[32];
  char gname[32];
  char devmajor[8];
  char devminor[8];
  char prefix[155];
  char padding[12];
};

union TarRecord {
  char buff[RECORDSIZE];
  TarHeader header;
};

size_t oct2size(char*);
void tgzReadFile(gzFile&, size_t, string&);

// read line of gnuzip file
int gzline(gzFile &, string &);

#endif // FILEOPT_H
