/*
 * Copyright (c) 2018  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2017-09-01 13:03:04
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-08 15:19:22
 */

#include "lineage.h"

/********************************************************************************
 * @brief Construct a new Lineage Handle:: Lineage Handle object.
 *        basic option of Lineage
 * @param taxfile 
 ********************************************************************************/
LineageHandle::LineageHandle(const string &taxfile) : taxfile(taxfile){
  setOutRankByTaxfile();
  setOutRank();
};

LineageHandle::LineageHandle(const string &taxfile, const string &revfile)
    : taxfile(taxfile), revfile(revfile) {
  setOutRankByTaxfile();
  setOutRank();
};

LineageHandle::LineageHandle(const string &taxfile, const string &revfile,
                             const string &abfile, const string &abtype)
    : taxfile(taxfile), revfile(revfile) {

  // revise the taxlevel if input give
  if (!abfile.empty()) {
    setRankByFile(abfile);
  } else {
    setOutRankByTaxfile();
  }

  // set out ranks
  if (!abtype.empty()) {
    outRankStr = abtype;
  }
  setOutRank();
};

void LineageHandle::setRankByFile(const string &file) {

  ifstream infile(file.c_str());
  if (!infile) {
    cerr << "Cannot found the input file " << file << endl;
    exit(5);
  }

  outRankStr.clear();
  for (string line; getline(infile, line);) {
    line = trim(line);
    if (line.empty()) {
      vector<string> words;
      separateWord(words, line);
      if (words.size() > 1) {
        rankmap[words[0]] = words[1].at(0);
        outRankStr += words[1].at(0);
      } else {
        rankmap[words[0]] = toupper(words[0].at(0));
        outRankStr += toupper(words[0].at(0));
      }
    }
  }
};


void LineageHandle::setOutRank() {

  map<char, string> abmap;
  for (auto &item : rankmap) {
    abmap[item.second] = item.first;
  };

  outrank.clear();
  for (const char &c : outRankStr) {
    if (abmap.find(c) != abmap.end()) {
      string str = abmap[c];
      str.front() = toupper(str.front());
      outrank.emplace_back(str, c);
    } else {
      cerr << "Can't find the full name of " << c << endl;
      exit(3);
    }
  }
};

size_t LineageHandle::nOutRanks() const { return outRankStr.size(); };

void LineageHandle::setOutRankByTaxfile() {
  ifstream tax(taxfile);
  if (tax) {
    outRankStr.clear();
    string line;
    do {
      getline(tax, line);
      trim(line);
    } while (line.empty());

    sregex_iterator iter(line.begin(), line.end(), prefix);
    sregex_iterator iterEnd;
    for (; iter != iterEnd; ++iter) {
      outRankStr += iter->str().at(1);
    }
    if (outRankStr.back() == 'T') {
      outRankStr.pop_back();
    }
    tax.close();
  }
};

/********************************************************************************
 * @brief obtained the lineage string from the taxfile
 * 
 * @param lngs a list of lineage for query and retrun
 ********************************************************************************/
void LineageHandle::getLineage(vector<Lineage> &lngs) {

  // use the cache file
  TaxMap taxmap;
  readTaxFile(taxmap);
  size_t nHit(0);
  for (auto &nm : lngs) {
    if (taxmap.find(nm.name) != taxmap.end()) {
      nm.name = taxmap[nm.name];
      nHit++;
    }
  }

  // output infomatics in taxfile
  theInfo("Find " + to_string(nHit) + "/" + to_string(lngs.size()) +
          " in file " + taxfile);

  // format the lineage and do revsion
  for (auto &lng : lngs) {
    format(lng.name);
  }

  // do revision
  Revision rev(revfile);
  if (!rev.empty()) {
    for (auto &lng : lngs) {
      rev.revise(lng.name);
    }
  }

  // check Repeat names
  checkRepeats(lngs);

  // set the unclassified items
  for (auto &nm : lngs) {
    if (nm.name.find(">" + undefStr + "<") != string::npos) {
      nm.def = false;
    }
  }
};

// read lineage string from input taxon file
void LineageHandle::readTaxFile(TaxMap &taxmap) {

  // get lineage map from file
  ifstream tax(taxfile);
  if (!tax) {
    cerr << "\nOpen " << taxfile
         << " failed, all lineage will obtained from sqlite" << endl;
  } else {

    // read the taxonomy map
    regex matchReg("([^; ]+)[; ]+([^; ]+)");
    for (string line; getline(tax, line);) {
      if (!line.empty()) {
        smatch matchs;
        regex_search(line, matchs, matchReg);
        if (!matchs.empty()) {
          string theName = matchs[1].str();
          string theTax = matchs[2].str();
          if (theTax.find("<T>") == string::npos)
            theTax += "<T>" + theName;
          taxmap[theName] = theTax;
        } else {
          cerr << "Error parse in " << line << " in file " << taxfile << endl;
        }
      }
    }
  }
  tax.close();
};

// format the lineage string according to outRank
void LineageHandle::format(string &atax) {

  string str = atax;
  atax.clear();
  for (const auto &t : outrank) {
    string tlab = "<";
    tlab += t.second;
    tlab += ">";
    auto bpos = str.find(tlab);
    if (bpos == string::npos) {
      atax += tlab;
      atax += undefStr;
    } else {
      auto epos = str.find_first_of('<', bpos + 1);
      atax += str.substr(bpos, epos - bpos);
    }
  }

  // for the strain name
  auto bpos = str.rfind("<T>");
  if (bpos == string::npos) {
    // don't find <T> label
    atax += "<T>";
    bpos = str.find_last_of('>');
    bpos = (bpos == string::npos) ? 0 : bpos + 1;
  }

  auto epos = str.find_first_of('<', bpos + 1);
  if (epos == string::npos) {
    atax += str.substr(bpos);
  } else {
    atax += str.substr(bpos, epos - bpos);
  }

  // for replace domain name to the missing kingdom name
  atax = regex_replace(atax, misKingdom, "<D>$1<K>$1");
};


// check the repeat items
void LineageHandle::checkRepeats(vector<Lineage> &lngs) {
  map<string, int> nameMap;
  for (auto &nm : lngs) {
    if (nameMap.find(nm.name) == nameMap.end()) {
      nameMap[nm.name] = 1;
    } else {
      theInfo("repeat strain: " + nm.name);
      nm.name += ".repeat-" + to_string(++nameMap[nm.name]);
    }
  }
}

/********************************************************************************
 * @brief convert the taxa to lineage string
 * 
 * @param lngs taxon level and taxon name
 * @return string a abbrivated string for the taxon
 ********************************************************************************/
string LineageHandle::lineageString(const vector<RankName> &lngs) const {
  string lineStr("");
  for (auto &lng : lngs) {
    lineStr += lineageString(lng);
  }
  return lineStr;
}

string LineageHandle::lineageString(const RankName &lng) const {
  string lineStr("");
  if (rankmap.find(lng.first) != rankmap.end()) {
    lineStr += "<";
    lineStr += rankmap.find(lng.first)->second;
    lineStr += ">";
    lineStr += lng.second;
  }
  return lineStr;
}

/********************************************************************************
 * @brief functions for the option on lineages
 *
 * @param taxstr
 * @param tax
 * @return size_t
 ********************************************************************************/
size_t parseLineage(const string &taxstr, vector<string> &tax) {
  size_t prev_pos = 0;
  size_t pos = 1;
  while ((pos = taxstr.find_first_of('<', pos)) != string::npos) {
    tax.emplace_back(taxstr.substr(0, pos));
    ++pos;
  }
  tax.emplace_back(taxstr);
  return tax.size();
};

size_t separateLineage(const string &taxstr, vector<string> &tax) {
  size_t prev_pos = 0;
  size_t pos = 1;
  while ((pos = taxstr.find_first_of('<', prev_pos + 1)) != string::npos) {
    tax.emplace_back(taxstr.substr(prev_pos, pos - prev_pos));
    prev_pos = pos;
  }
  tax.emplace_back(taxstr.substr(prev_pos));
  return tax.size();
};

string lastName(const string &str) {
  return str.substr(str.find_last_of('<'));
};

string lastNameNoRank(const string &str) {
  size_t pos = str.find_last_of('>');
  if (pos == string::npos)
    return str;
  return str.substr(pos + 1);
};

string commonLineage(const string &lngA, const string &lngB) {

  string comlng;
  // check whether at root/empty
  if (lngA.empty() || lngB.empty()) {
    return comlng;
  }

  // compare lineages
  vector<string> lngListA;
  size_t nlngA = parseLineage(lngA, lngListA);
  vector<string> lngListB;
  size_t nlngB = parseLineage(lngB, lngListB);

  size_t nlev = nlngA < nlngB ? nlngA : nlngB;
  for (size_t i = 0; i < nlev; ++i) {
    if (lngListA[i] == lngListB[i]) {
      comlng = lngListA[i];
    } else {
      break;
    }
  }
  return comlng;
};

// common lineage for a list
string commonLineage(const vector<string> &lngs) {
  string comlng;
  if (!lngs.empty()) {
    auto iter = lngs.begin();
    comlng = *iter;
    for (++iter; iter != lngs.end(); ++iter) {
      comlng = commonLineage(comlng, *iter);
    }
  }
  return comlng;
};

// number ranks
size_t nRanks(const string &lng) {
  size_t n = 0;
  for (auto c : lng) {
    if (c == '<')
      n++;
  }
  return n;
}
