#ifndef FUNCTIONS_HPP_
#define FUNCTIONS_HPP_
#include "tools.hpp"
#include "TTree.h"
#include <map>
#include <vector>
#include <iostream>
#include <string>

using namespace std;

bool isFake(TTree* event);
bool doTTbarCombination(TTree* event);
string doHeavyFlavor(TTree* event);
std::map<string, float> GetSysWeights(TTree* event);

#endif
