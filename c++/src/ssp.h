#include <vector>
#include <map>
#include <Sequence/SimpleSNP.hpp>
using namespace std;

typedef vector<uint32_t> IntVec;
typedef vector<string> StrVec;
typedef map< uint32_t, string > IntStrMap;
typedef map< string, string > StrStrMap;
typedef map< string, vector<uint32_t> > StrVecMap;

StrVec get_acc_ids(const string& f_id, const string& opt); 
IntVec get_idx(StrVec& si, StrVec& sa); 
IntVec sample_serial(const uint32_t& n, uint32_t& m); 
IntVec sample_random(const uint32_t& n, const uint32_t& m); 

namespace Sequence {
    class SSP: public SimpleSNP {
      private:
      public:
        unsigned nind, npos;
        std::vector<std::string> labels;

        std::istream& read (std::istream& s) throw (Sequence::badFormat, std::exception);
        std::ostream& write (std::ostream& o, const string& fmt, const string& chr="chr0"); 
        std::ostream& write_expand(std::ostream& o, const string& dir, const string& seqid);
        void set_label(unsigned i, const string& label);
        unsigned get_label_idx(const string& label);

        void ApplyFreqFilter(unsigned mincount, bool haveOutgroup=false, unsigned outgroup=0);
        void RemoveMultiHits(bool skipOutgroup=false, unsigned outgroup=0);
        void RemoveMissing(bool skipOutgroup=false, unsigned outgroup=0);
        void RemoveAmbiguous(bool skipOutgroup=false, unsigned outgroup=0);

        void RemoveMono(bool skipOutgroup=false, unsigned outgroup=0);
        void FilterMissing(const int& co_mis, bool skipOutgroup=false, unsigned outgroup=0);
        void subset(IntVec& idxs_ind, IntVec& idxs_pos);
    };
}


