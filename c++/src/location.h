#include <string>
#include <vector>
using namespace std;

typedef vector<int> IntVec;
typedef set<int> IntSet;
struct Location {
    uint32_t beg;
    uint32_t end;
    double score;
    int idx_ext;
    IntSet idxs;
};
typedef vector< Location > LocVec;

bool compare_loc (const Location& l1, const Location& l2);
LocVec tiling(const LocVec& lvi,  const bool& flag_max = false);

