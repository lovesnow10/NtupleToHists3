/*
  A simple DataSets Handler, support some basic funtions
*/
#ifndef DSHANDLER_HPP_
#define DSHANDLER_HPP_

#include <iostream>
#include <string>
#include <map>
#include <vector>

// Root Headers
#include "TSystemDirectory.h"
#include "tools.hpp"

using namespace std;

class DSHandler {
private:
  std::vector<string> mPathsVec;
  std::map<int, std::vector<string>> mFilesMap;
  std::map<int, string> mSampleTypeMap;

  int mTotalDS;
  int mCurrentDS;

  std::vector<string> GrabRootFiles(string path);

  std::vector<string> tmpNone;

  string GetSampleID(string path);

public:
  DSHandler(){};
  DSHandler(string InputTxt);
  virtual ~DSHandler(){};

  const std::vector<string> &GetFiles(int nDS) const;

  void AddPath(string path);
  void Clear();
  void DeletePath(string path);
  bool Initialize();
  void PrintCurrentDSInfo() const
  {
    PrintDSInfo(mCurrentDS);
  };
  void PrintDSInfo(int nDS) const;

  const std::vector<string> &Next();
  string GetSampleType(int nDS) {return mSampleTypeMap.at(nDS);};
  int GetSampleIndex() {return this->mCurrentDS;};
  int GetNSamples() {return mTotalDS;};
};

#endif
