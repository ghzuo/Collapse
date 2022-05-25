/*
 * Copyright (c) 2022  Wenzhou Institute, University of Chinese Academy of
 * Sciences. See the accompanying Manual for the contributors and the way to
 * cite this work. Comments and suggestions welcome. Please contact Dr.
 * Guanghong Zuo <ghzuo@ucas.ac.cn>
 *
 * @Author: Dr. Guanghong Zuo
 * @Date: 2022-03-16 12:10:27
 * @Last Modified By: Dr. Guanghong Zuo
 * @Last Modified Time: 2022-05-25 01:57:46
 */

#include "lineage.h"

ostream &operator<<(ostream &os, const Lineage &lng) {
  vector<string> ranks;
  separateLineage(lng.name, ranks);
  for (auto &str : ranks) {
    str = lastNameNoRank(str);
  }
  os << "[\"";
  os << ranks.back() << "\",\"";
  os << strjoin(ranks.begin(), ranks.end() - 1, "\",\"");
  os << "\"]";
  return os;
};

/********************************************************************************
 * @brief Construct a new Lineage Data:: Lineage Map.
 *        basic option of Lineage
 * @param taxfile
 ********************************************************************************/
LngData::LngData() { rank = rank->create(); };

LngData::LngData(const string &dbpath, const string &taxfile,
                 const string &revfile)
    : tdbpath(dbpath), revfile(revfile) {
  separateWord(tflist, taxfile);
  rank = rank->create();
};

/********************************************************************************
 * @brief obtained the lineage string from the taxfile
 *
 * @param lngs a list of lineage for query and retrun
 ********************************************************************************/
void LngData::getLineage(vector<string> &nmlist) {

  // initial the data
  for (auto &nm : nmlist) {
    data.emplace_back(nm);
  }

  // find the lineage in the taxfiles
  size_t nHit = getLngFromFile();

  // find lineage in database
  if (nHit < data.size() && fileExists(tdbpath))
    getLngFromDB();

  // format the lineage and do revsion
  for (auto &lng : data) {
    rank->format(lng.name);
  }

  // do revision
  if (!revfile.empty()) {
    Revision rev(revfile);
    if (!rev.empty()) {
      int nReplace(0);
      for (auto &lng : data) {
        if (rev.revise(lng.name) > 0)
          ++nReplace;
      }
      theInfo("There are " + to_string(nReplace) +
              " lineage had been changed by revision file");
    }
  }

  // check Repeat names
  checkRepeats();

  // set the unclassified items
  for (auto &lng : data) {
    lng.def = rank->wellDefined(lng.name);
  }
};

// read the lineage string from input taxon file
size_t LngData::getLngFromFile() {

  size_t nTotal(0);
  for (auto &taxfile : tflist) {
    ifstream tax(taxfile);
    if (!tax) {
      theInfo("Open file " + taxfile + " failed, skip this file for lineage");
    } else {
      TaxMap taxmap;
      string format = toLower(getsuffix(taxfile));
      if (format.compare("lns") == 0) {
        readListFile(tax, taxmap);
      } else if (taxfile.compare("json") == 0) {
        readJsonFile(tax, taxmap);
      } else {
        readCsvFile(tax, taxmap);
      }

      size_t nHit(0);
      if (!taxmap.empty()) {
        for (auto &nm : data) {
          if (!nm.def) {
            auto iter = taxmap.find(nm.name);
            if (iter != taxmap.end()) {
              nm.name = rank->lineage(iter->second, nm.name);
              nHit++;
              nm.def = true;
            }
          }
        }
      }
      theInfo("Found " + to_string(nHit) + "/" + to_string(data.size()) +
              " lineages in file " + taxfile);
      nTotal += nHit;
    }
    tax.close();
  }

  return nTotal;
};

// get lineage from database file
void LngData::getLngFromDB() {
  TaxaDB taxadb(tdbpath);
  size_t mHit(0);
  for (auto &lng : data) {
    if (!lng.def) {
      string res = taxadb.search(lng.name);
      if (!res.empty()) {
        lng.name = rank->lineage(res, lng.name);
        mHit++;
      }
    }
  }
  theInfo("Found " + to_string(mHit) + "/" + to_string(data.size()) +
          " lineages in taxa database " + tdbpath);
}

// read lineage string from list file
void LngData::readListFile(istream &tax, TaxMap &taxmap) {
  // get the content of file
  vector<string> lines;
  for (string line; getline(tax, line);) {
    if (!line.empty()) {
      lines.emplace_back(line);
    }
  }

  // reset the out ranks
  string outRankStr;
  regex prefix{rank->addMark("[A-Za-z]")};
  sregex_iterator iter(lines.front().begin(), lines.front().end(), prefix);
  sregex_iterator iterEnd;
  for (; iter != iterEnd; ++iter) {
    outRankStr += iter->str().at(1);
  }
  if (outRankStr.back() == 'T') {
    outRankStr.pop_back();
  }
  rank->setOutRank(outRankStr);

  // read the taxonomy map
  for (auto &line : lines) {
    if (!line.empty()) {
      vector<string> items;
      separateWord(items, line, " ,;");
      if (items.size() > 1) {
        taxmap[items[0]] = items[1];
      } else {
        string theName = lastName(items[0]);
        if (theName.at(1) == 'T') {
          taxmap[lastNameNoRank(items[0])] = items[0];
        } else {
          theInfo("Error parse in " + line);
        }
      }
    }
  }
};

// read lineage string from csv file
void LngData::readCsvFile(istream &tax, TaxMap &taxmap) {
  // get the content of file
  vector<vector<string>> content;
  for (string line; getline(tax, line);) {
    if (!line.empty()) {
      vector<string> words;
      separateWord(words, line);
      content.emplace_back(words);
    }
  }

  // set out rank and symbol by header
  vector<char> symlist;
  vector<pair<string, char>> theRanks;
  vector<string> items(content.front().begin() + 1, content.front().end());
  for (auto &item : items) {
    char symbol;
    string level;
    size_t pos = item.find_first_of("(:");
    if (pos == string::npos) {
      level = item;
      symbol = rank->getSymbol(level);
      if (symbol == ' ') {
        cerr << "Cannot find the symbol of " << level << endl;
        exit(10);
      }
    } else {
      symbol = item.at(pos + 1);
      level = item.substr(0, pos);
    }
    symlist.emplace_back(symbol);
    theRanks.emplace_back(level, symbol);
  }
  rank->setOutRank(theRanks);

  // get lineage data
  for (auto iter = content.begin() + 1; iter != content.end(); ++iter) {
    string name = iter->front();
    string taxstr;
    for (size_t i = 0; i < symlist.size(); ++i) {
      taxstr += TaxaRank::addMark(symlist[i]);
      taxstr += iter->at(i + 1);
    }
    taxmap[name] = taxstr;
  }
};

// read lineage string from json file
void LngData::readJsonFile(istream &tax, TaxMap &taxmap) {
  Json taxjson = Json::parse(tax);

  // get the symbol list
  vector<char> symlist;
  vector<pair<string, char>> theRanks;
  for (auto &lvl : taxjson.at("rank")) {
    char symbol = lvl.at("symbol").get<std::string>().front();
    symlist.emplace_back(symbol);
    theRanks.emplace_back(lvl.at("level"), symbol);
  }
  rank->setOutRank(theRanks);

  // get the lineage records
  for (auto &jlng : taxjson.at("lineage")) {
    string name = jlng.at(0);
    string taxstr;
    for (size_t i = 0; i < symlist.size(); ++i) {
      taxstr += TaxaRank::addMark(symlist[i]);
      taxstr += jlng.at(i + 1);
    }
    taxmap[name] = taxstr;
  }
};

// check the repeat items
void LngData::checkRepeats() {
  map<string, int> nameMap;
  for (auto &nm : data) {
    if (nameMap.find(nm.name) == nameMap.end()) {
      nameMap[nm.name] = 1;
    } else {
      theInfo("Found repeat strain: " + nm.name);
      nm.name += ".repeat-" + to_string(++nameMap[nm.name]);
    }
  }
}

Lineage &LngData::operator[](size_t ndx) { return data[ndx]; };
size_t LngData::size() const { return data.size(); };
void LngData::getStat(vector<pair<int, int>> &stat) const {
  // initial the stat
  stat.resize(rank->nOutRanks());
  for (auto &item : stat) {
    item = make_pair(0, 0);
  }

  // do the statistics
  vector<set<string>> uniqTaxa(rank->nOutRanks());
  for (auto &lng : data) {
    vector<string> theLng;
    separateLineage(lng.name, theLng);
    for (size_t i = 0; i < rank->nOutRanks(); ++i) {
      if (rank->wellDefined(theLng[i])) {
        uniqTaxa[i].insert(theLng[i]);
      } else {
        ++stat[i].second;
      }
    }
  }

  // the the number unique taxon at every level
  for (size_t i = 0; i < rank->nOutRanks(); ++i) {
    stat[i].first = uniqTaxa[i].size();
  }
};

void LngData::output(const string &fname, const string &format) const {
  ofstream LNG(fname.c_str());
  if (!LNG) {
    cerr << "Open file failed at " << fname << endl;
    exit(10);
  }

  string outformat = toLower(format);
  if (outformat.compare("json") == 0) {
    stringstream jsonbuf;
    outjson(jsonbuf);
    theJson(LNG, jsonbuf.str());
    LNG << endl;
  } else if (outformat.compare("lns") == 0) {
    outlist(LNG);
  } else {
    outcsv(LNG);
  }
};

void LngData::outjson(ostream &os) const {
  // output lineage of leaf;
  os << "{\"lineage\":";
  os << "[" << strjoin(data.begin(), data.end(), ",") << "]";
  os << ",";

  // get the statistics of lineage
  vector<pair<int, int>> stat;
  getStat(stat);

  // output the ranks and statistics
  vector<string> jslist;
  for (int i = 0; i < rank->nOutRanks(); ++i) {
    jslist.emplace_back("{\"level\":\"" + rank->outrank[i].first +
                        "\",\"symbol\":\"" + rank->outrank[i].second +
                        "\",\"nTaxa\":\"" + to_string(stat[i].first) +
                        "\",\"nUndef\":\"" + to_string(stat[i].second) + "\"}");
  }
  os << "\"rank\":[" << strjoin(jslist.begin(), jslist.end(), ',') << "]}";
};

void LngData::outcsv(ostream &os) const {
  // the table head
  os << "Strain(T),";
  rank->outRanksCSV(os);

  for (auto &lng : data) {
    vector<string> ranks;
    separateLineage(lng.name, ranks);
    for (auto &str : ranks) {
      str = lastNameNoRank(str);
    }
    os << "\n" << ranks.back() << ",";
    os << strjoin(ranks.begin(), ranks.end() - 1, ",");
  }
};

void LngData::outlist(ostream &os) const {
  vector<string> lines;
  for (auto &lng : data) {
    lines.emplace_back(lastNameNoRank(lng.name) + " " + lng.name);
  }
  os << strjoin(lines.begin(), lines.end(), "\n");
};
