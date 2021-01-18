/*
 * Copyright (c) 2020  T-Life Research Center, Fudan University, Shanghai,
 * China. See the accompanying Manual for the contributors and the way to cite
 * this work. Comments and suggestions welcome. Please contact Dr. Guanghong Zuo
 * <ghzuo@fudan.edu.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2020-12-09 16:10:07
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2020-12-24 15:22:26
 */

#include "taxarank.h"

/********************************************************************************
 * @brief factory function for TaxaRank to keep it only one copy in program
 *
 ********************************************************************************/
pair<char, char> TaxaRank::mark = make_pair('<', '>');

TaxaRank *TaxaRank::pRank = nullptr;
TaxaRank *TaxaRank::create() {
  if (pRank == nullptr) {
    pRank = new TaxaRank();
  }
  return pRank;
};

/********************************************************************************
 * @brief set rank
 *
 ********************************************************************************/
void TaxaRank::initial(const string &taxfile, const string &abfile,
                       const string &abtype) {
  // revise the taxlevel if input give
  string outRankStr = setRankByFile(abfile);

  // set out ranks
  if (!abtype.empty()) {
    outRankStr = abtype;
  } else if (outRankStr.empty()) {
    outRankStr = setOutRankByTaxfile(taxfile);
  }

  setOutRank(outRankStr);
};

void TaxaRank::initial(const string &abfile, const string &abtype) {
  // revise the taxlevel if input give
  string outRankStr = setRankByFile(abfile);

  // set out ranks
  if (!abtype.empty()) {
    outRankStr = abtype;
  }
  setOutRank(outRankStr);
};

string TaxaRank::setRankByFile(const string &file) {

  string outRankStr;
  if (!file.empty()) {

    ifstream infile(file.c_str());
    if (!infile) {
      cerr << "Cannot found the input file " << file << endl;
      exit(5);
    }

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
  }
  return outRankStr;
};

string TaxaRank::setOutRankByTaxfile(const string &file) {
  string outRankStr;
  ifstream tax(file);
  if (tax) {
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
  return outRankStr;
};

void TaxaRank::setOutRank(const string &outRankStr) {

  if (!outRankStr.empty()) {
    map<char, string> abmap;
    for (auto &item : rankmap) {
      abmap[item.second] = item.first;
    };
    abmap['D'] = "Domain";

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
  }
};

void TaxaRank::outRanksJson(ostream &os) {
  vector<string> jslist;
  for (auto &rank : outrank) {
    jslist.emplace_back("{\"level\":\"" + rank.first + "\",\"symbol\":\"" +
                        rank.second + "\"}");
  }
  os << "[" << strjoin(jslist.begin(), jslist.end(), ',') << "]";
};

/********************************************************************************
 * @brief for check ranks
 *
 ********************************************************************************/
size_t TaxaRank::nOutRanks() const { return outrank.size(); };

bool TaxaRank::wellDefined(const string &lng) {
  return lng.find(undefSym) == string::npos;
};

string TaxaRank::lineage(const string &lng, const string &str) {
  if (lng.find(strainMark) == string::npos) {
    return lng + strainMark + str;
  }
  return lng;
};

// format the lineage string according to outRank
void TaxaRank::format(string &atax) {

  string str = atax;
  atax.clear();
  for (auto &t : outrank) {
    string tlab = addMark(t.second);
    auto bpos = str.find(tlab);
    if (bpos == string::npos) {
      atax += tlab;
      atax += undefStr;
    } else {
      auto epos = str.find_first_of(mark.first, bpos + 1);
      atax += str.substr(bpos, epos - bpos);
    }
  }

  // for the strain name
  auto bpos = str.rfind(strainMark);
  if (bpos == string::npos) {
    // don't find <T> label
    atax += strainMark;
    bpos = str.find_last_of(mark.second);
    bpos = (bpos == string::npos) ? 0 : bpos + 1;
  }

  auto epos = str.find_first_of(mark.first, bpos + 1);
  if (epos == string::npos) {
    atax += str.substr(bpos);
  } else {
    atax += str.substr(bpos, epos - bpos);
  }

  // for replace domain name to the missing kingdom name
  atax = regex_replace(atax, misKingdom, comKingdom);
};

/********************************************************************************
 * @brief convert the taxa to lineage string
 *
 * @param lngs taxon level and taxon name
 * @return string a abbrivated string for the taxon
 ********************************************************************************/
string TaxaRank::rankString(const vector<RankName> &rnlist) const {
  string lineStr("");
  for (auto &rn : rnlist) {
    lineStr += rankString(rnlist);
  }
  return lineStr;
}

string TaxaRank::rankString(const RankName &rn) const {
  string lineStr("");
  if (rankmap.find(rn.first) != rankmap.end()) {
    lineStr += mark.first;
    lineStr += rankmap.find(rn.first)->second;
    lineStr += mark.second;
    lineStr += rn.second;
  }
  return lineStr;
}

/********************************************************************************
 * @brief functions for the option on lineages
 *
 ********************************************************************************/
size_t parseLineage(const string &taxstr, vector<string> &tax) {
  size_t prev_pos = 0;
  size_t pos = 1;
  while ((pos = taxstr.find_first_of(TaxaRank::mark.first, pos)) !=
         string::npos) {
    tax.emplace_back(taxstr.substr(0, pos));
    ++pos;
  }
  tax.emplace_back(taxstr);
  return tax.size();
};

size_t separateLineage(const string &taxstr, vector<string> &tax) {
  size_t prev_pos = 0;
  size_t pos = 1;
  while ((pos = taxstr.find_first_of(TaxaRank::mark.first, prev_pos + 1)) !=
         string::npos) {
    tax.emplace_back(taxstr.substr(prev_pos, pos - prev_pos));
    prev_pos = pos;
  }
  tax.emplace_back(taxstr.substr(prev_pos));
  return tax.size();
};

string lastName(const string &str) {
  return str.substr(str.find_last_of(TaxaRank::mark.first));
};

string lastNameNoRank(const string &str) {
  size_t pos = str.find_last_of(TaxaRank::mark.second);
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

// delete the strain section in lineage
string delStrain(const string &lng) {
  size_t npos = lng.find(TaxaRank::addMark('T'));
  if (npos != string::npos) {
    return lng.substr(0, npos);
  }
  return lng;
};

// number ranks
size_t nRanks(const string &lng) {
  size_t n = 0;
  for (auto c : lng) {
    if (c == TaxaRank::mark.first)
      n++;
  }
  return n;
}
