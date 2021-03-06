/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2018-08-01 22:28:21
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-08 14:59:37
 */

#define THEINFO
#include "info.h"
Info theInfo;

/********************************************************************************
 * @brief member for Timer
 * 
 ********************************************************************************/
Timer::Timer() : start(std::chrono::system_clock::now()){};

double Timer::elapsed() {
  auto now = std::chrono::system_clock::now();
  double time =
      std::chrono::duration_cast<std::chrono::milliseconds>(now - start)
          .count();
  return time / 1000.0;
};

/********************************************************************************
 * @brief Construct a new Info:: Info object
 * 
 ********************************************************************************/
Info::Info() : dep(0), quiet(false){};

Info::~Info() {
  //... destroy the lock
  if (!quiet)
    cerr << "*** ALL Section: Complete Program, Time Elapsed: "
         << mytimer.elapsed() << " s" << endl;
}

void Info::operator()(const string &str, int idep) {
  if (!quiet) {
    if (idep > 0)
      dep += idep;
    string indent(dep * 4, ' ');
    if (idep < 0)
      dep += idep;
    indent += "*** ";

    string sstr(str);
    size_t npos(0);
    while ((npos = sstr.find('\n', npos)) != std::string::npos) {
      ++npos;
      sstr.insert(npos, indent);
    }
    cerr << indent << sstr << ", Time Elapsed: " << mytimer.elapsed() << " s"
         << endl;
  }
};