#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>
using namespace std;

typedef pair< uint32_t, uint32_t > LocPair;
typedef vector< LocPair > LocPairVec;
typedef vector<int> IntVec;
typedef set<int> IntSet;
typedef vector< IntSet > IntSetVec;
typedef vector<double> DubVec;

void pos_split(const LocPair& lpi, const DubVec& scores, LocPair& lpo, IntVec& idxs, const int& opt = 1);

#endif
