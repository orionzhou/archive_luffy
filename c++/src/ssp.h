#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>
#include <Sequence/PolyTable.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/SeqExceptions.hpp>
using namespace std;
using boost::format;
using namespace boost::assign;

typedef vector<uint32_t> IntVec;
typedef vector<string> StrVec;
typedef map< uint32_t, string > IntStrMap;
typedef map< string, string > StrStrMap;
typedef map< string, vector<uint32_t> > StrVecMap;

StrVec get_acc_ids(const string& f_id, const string& opt); 
IntVec get_acc_idx(const string& f_id, const string& opt_s, StrVec& ids, const int& opt); 
IntVec sample_serial(const uint32_t& n, uint32_t& m); 
IntVec sample_random(const uint32_t& n, const uint32_t& m); 

namespace Sequence {
  class SimpleSNP:public PolyTable {
    private:
      mutable bool Diploid,isoFemale;
      bool haveOutgroup;
      std::vector<std::string> _names;
    public:
      SimpleSNP (const bool diploid =0,const bool isofemale=0)  : PolyTable(),
          Diploid(diploid),isoFemale(isofemale),haveOutgroup(false)
      {}
      ~ SimpleSNP (void) {}
     
      bool outgroup(void) const;
      std::string label(unsigned i) const;
      void set_label(unsigned i, const string & str);
      
      std::istream & read (std::istream & s) throw (Sequence::badFormat,std::exception);
      std::ostream & print(std::ostream & o) const;
      void write(const string & fo, const string & out_format, const string & chr="chr0"); 
      void write_expand(const string & dirO, const string& seqid);
      
      void RemoveMono(bool skipOutgroup=false, unsigned outgroup=0);
      void FilterMissing(const int & co_mis, bool skipOutgroup=false, unsigned outgroup=0);
      void subset(IntVec & idxs_ind, IntVec & idxs_pos);
      void set_outgroup(const string & ind_og);
      void clean_outgroup(unsigned & n_N, unsigned & n_fix);
    };
    typedef SimpleSNP Hudson2001;
}


