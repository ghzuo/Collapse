/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-12-12 10:07:59
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-12 13:43:57
 */

#include "fileOpt.h"

/********************************************************************************
 * @brief Options on tar and zlib files
 *
 ********************************************************************************/

size_t oct2size(char *cstr) {
  size_t nsize(0);
  size_t nbase = 1;
  for (int i = 10; i > 0; i--) {
    nsize += nbase * (cstr[i] - 48);
    nbase *= 8;
  }

  return nsize;
};

void tgzReadFile(gzFile &fp, size_t nsize, string &str) {
  char buf[RECORDSIZE];
  for (size_t i = 0; i < nsize; i += RECORDSIZE) {
    gzread(fp, buf, sizeof(buf));
    str.append(buf, sizeof(buf));
  }
  str.resize(nsize);
  str.shrink_to_fit();
};

// get line from gz file
int gzline(gzFile &fp, string &line) {
  line.clear();
  char ch = gzgetc(fp);
  for (; ch != -1; ch = gzgetc(fp)) {
    if (ch == '\n') {
      return ch;
    } else {
      line += ch;
    }
  }
  return ch;
};