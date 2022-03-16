/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of Sciences.
 * See the accompanying Manual for the contributors and the way to cite this work.
 * Comments and suggestions welcome. Please contact
 * Dr. Guanghong Zuo <ghzuo@ucas.ac.cn>
 * 
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-03-16 12:28:26
 */

#include "reviseList.h"

Revision::Revision(const string &file) {

  /// get the revsion items
  ifstream is(file);
  if (is) {
    for (string line; getline(is, line);) {
      if (!line.empty()) {
        string str(trim(line.substr(0, line.find_first_of('#'))));
        if (!str.empty()) {
          vector<string> w;
          if (separateWord(w, str) == 2) {
            chglist.emplace_back(w[0], w[1]);
          } else {
            vector<size_t> vpos;
            size_t prev_pos = 0;
            size_t pos = 0;
            while ((pos = str.find_first_of('/', pos)) != string::npos)
              vpos.emplace_back(++pos);

            if (vpos.size() == 3) {
              chglist.emplace_back(str.substr(vpos[0], vpos[1] - vpos[0] - 1),
                                   str.substr(vpos[1], vpos[2] - vpos[1] - 1));
            } else {
              cerr << "Error parse in " << str << endl;
              exit(2);
            }
          }
        }
      }
    }
    is.close();

    if (chglist.empty()) {
      theInfo("There are no revsion on lineage in file " + file);
    } else {
      theInfo("There are " + to_string(chglist.size()) +
              " revsions on lineages in file " + file);
    }
  } else {
    theInfo("Open file " + file + ", no revsion on lineage");
  }
}

bool Revision::empty() const { return chglist.empty(); };

void Revision::revise(vector<string> &nmlist) {
  for (str2str &p : chglist)
    for (auto &nm : nmlist)
      nm = regex_replace(nm, regex(p.first), p.second);
}

void Revision::revise(string &nm) {
  for (str2str &p : chglist)
    nm = regex_replace(nm, regex(p.first), p.second);
}