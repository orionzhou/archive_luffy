#include <Sequence/PolyTable.hpp>
#include <Sequence/stateCounter.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/Portability/StringStreams.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <cstring>
#include <cctype>
#include <cassert>
#include "ssp.h" 
using namespace std;
using namespace boost::assign;
using std::vector;
using std::string;
using boost::format;

template< class T >
struct next {
    T seed;
    next( T seed ) : seed(seed) { }
    T operator()() {
        return seed++;
    }
};
StrVec get_acc_ids(const string& f_id, const string& opt) {
    ifstream fhi(f_id.c_str());
    if(!fhi.is_open()) {
        cout << format("cannot read file: %s\n") % f_id; exit(1);
    }
    int idx = -1;
    string line;

    getline(fhi, line);
    vector<string> ss;
    boost::split(ss, line, boost::is_any_of("\t"));
    for(uint32_t i = 0; i < ss.size(); i ++)
        if(ss[i] == opt) idx = i;
    if(idx == -1) {
        cout << "unknown opt: " << opt << endl;
        exit(1);
    }

    StrVec ids;
    while(fhi.good()) {
        getline(fhi, line);
        boost::split(ss, line, boost::is_any_of("\t"));
        if(ss[idx] == "1") ids.push_back(ss[0]);
    }
    return ids;
}
IntVec get_acc_idx(const string& f_id, const string& opt_s, StrVec& ids_s, const int& opt) {
    StrVec ids;
    if( opt == 1 ) {
        ids = get_acc_ids(f_id, "acc85"); 
    } else if ( opt == 2 ) {
        ids = get_acc_ids(f_id, "acc289"); 
    } else {
        cout << "unknown ind opt: " << opt << endl;
        exit(1);
    }
    ids_s = get_acc_ids(f_id, opt_s);
    cout << format("  selecting ids: %s\n") % boost::join(ids_s, " ");

    map<string, int> tmp;
    for(uint32_t i = 0; i < ids.size(); i ++) tmp[ids[i]] = i;

    IntVec idxs;
    for(uint32_t i = 0; i < ids_s.size(); i ++) {
        string id = ids_s[i];
        map<string, int>::iterator it = tmp.find(id);
        if(it == tmp.end()) {
            cout << format("  ind[%s] not found in opt[%d]\n") % id % opt;
            exit(1);
        }
        idxs.push_back( it->second );
        //    cout << format("%s: %d\n") % id % it->second;
    }
    return idxs;
}
IntVec sample_serial(const uint32_t& n, uint32_t& m) {
    IntVec iv;
    if(m > n || m <= 0) m = n;
    double increment = n / m;
    for(uint32_t i = 0; i < m; i ++) 
        iv.push_back(1 + int(i * increment) - 1);
    return iv;
}
IntVec sample_random(const uint32_t& n, const uint32_t& m) {
    IntVec iv;
    for(uint32_t i = 0; i < m; i ++) 
        iv.push_back(1 + int(n * rand()/(RAND_MAX+1.0)) - 1);
    return iv;
}

namespace Sequence {
  std::istream & SimpleSNP::read (std::istream & s)  throw (Sequence::badFormat,std::exception)
  {
    unsigned nsam, nsites, i, j;
    char ch;
    if(!(s >> nsam >> nsites))
      throw badFormat("SimpleSNP.cc: file did not start with s nsam nsites");
    if (Diploid)
      nsam *= 2;
    std::vector<double> _positions(nsites);

    for (i = 0; i < nsites; ++i)
      {
        if (! (s >> _positions[i] >> std::ws))
          throw badFormat("SimpleSNP.cc: error in processing site positions");
      }
    std::string outgroup,temp,temp2;
    std::getline(s,temp);
    istr check(temp),check2(temp);
    unsigned nc = 0;
    while (! check.eof() )
      {
	check >> temp2 >> std::ws;
	++nc;
      }
    bool haveOGlabel = (nc == nsites + 1) ? true : false;
    _names.resize(nsam+1);
    if(haveOGlabel)
      check2 >> _names[0];
    else
      _names[0]="anc";
    outgroup.resize(nsites);
    for (i = 0; i < nsites; ++i)
      {
 	if(!(check2 >> ch))
 	  {
 	    throw badFormat("SimpleSNP.cc: error reading in seg. sites");
 	  }
        ch = char(toupper(ch));
        switch (ch)
          {
          case '?':
            outgroup[i] = 'N';
            break;
          default:
            outgroup[i] = ch;
            break;
          }
      }
    outgroup[i] = '\0';
    unsigned numn = 0;
    //count # of ambiguous bases in outgroup
    for (i = 0; i <  outgroup.length(); ++i)
      {
        if (char(std::toupper(outgroup[i])) == 'N')
          ++numn;
      }
    unsigned haveOG = 0;

    std::vector<std::string> _data;
    if(numn == outgroup.length())
      //if all outgroup characters are ambiguous,
      //then we have no outgroup, so we don't include
      //it in the data vector, and so we don't allocate space
      {
        _data.resize(nsam);
        for (i = 0 ; i < nsam ; ++i)
          _data[i].resize(nsites);
      }
    else
      //otherwise,we have an outgroup and need to allocate space
      {
        _data.resize(nsam+1);
        for (i = 0 ; i < nsam+1 ; ++i)
          _data[i].resize(nsites);

        _data[0] = outgroup;
        haveOG=1;
        haveOutgroup = true;
      }

    for (i = 0 + haveOG; i < nsam + haveOG; ++i)
      {
        string name;	//don't store the name
        if(!(s >> name))
          throw badFormat("SimpleSNP.cc: error processing sequences");
        //_names[i-(haveOG+haveOGlabel)] = name;
	_names[i-haveOG+1] = name;
        char *temp = new char[nsites+1];
        char *temp2 = NULL;
        if (Diploid)
          {
            //_names[i-(haveOG+haveOGlabel)+1] = name;
	    _names[i-haveOG+2] = name;
            temp2 = new char[nsites+1];
          }
        for (j = 0; j < nsites; ++j)
          {
            if(!(s >> ch))
              throw badFormat("SimpleSNP.cc: error processing sequenes");
            ch = char(toupper(ch));
            if (Diploid)
              {
                switch (char(std::toupper(ch)))	//use IUPAC ambiguity symbols
                  {
                  case '?':
                    temp[j] = 'N';
                    temp2[j] = 'N';
                    break;
                  case 'M':
                    temp[j] = 'A';
                    temp2[j] = 'C';
                    break;
                  case 'R':
                    temp[j] = 'A';
                    temp2[j] = 'G';
                    break;
                  case 'W':
                    temp[j] = 'A';
                    temp2[j] = 'T';
                    break;
                  case 'S':
                    temp[j] = 'C';
                    temp2[j] = 'G';
                    break;
                  case 'Y':
                    temp[j] = 'C';
                    temp2[j] = 'T';
                    break;
                  case 'K':
                    temp[j] = 'G';
                    temp2[j] = 'T';
                    break;
                  default:
                    temp[j] = ch;
                    temp2[j] = ch;
                  }
              }
            else if (isoFemale)
              {
                //not implemented yet!!
              }
            else
              {
                switch (ch)
                  {
                  case '?':
                    temp[j] = 'N';
                    break;
                  default:
                    temp[j] = ch;
                    break;
                  }
              }
          }
        temp[j] = '\0';
        if (Diploid)
          {
            temp2[j] = '\0';
            _data[i++] = temp;
            _data[i] =  temp2;
            delete [] temp2;
          }
        else
          _data[i] = temp;

        delete [] temp;
      }
    //assign the data to the base class
    if (_data.size() != nsam+haveOG)
      {
	throw (Sequence::badFormat("SimpleSNP::read() -- number of sequences does not match input value"));
      }     
    PolyTable::assign(&_positions[0],_positions.size(),&_data[0],_data.size());
    RemoveInvariantColumns(this);
    return s;
  }
  std::ostream & SimpleSNP::print(std::ostream &o) const
  {
    o << this->size()-this->haveOutgroup <<'\t' <<this->numsites() << '\n';

    for(unsigned i = 0 ; i < this->numsites()-1 ; ++i) {
    	o << int(this->position(i)) << '\t';
    }
    o << int(this->position(this->numsites()-1)) << '\n';

    if (haveOutgroup == true) {
    	if(_names.empty()) {
  	    o << "anc";
      } else {
  	    o << _names[0];
	    }
      for(unsigned i = 0 ; i < this->numsites() ; ++i) {
        o << '\t' << (*this)[0][i];
      }
      o << '\n';
      for(unsigned i = 1 ; i < this->size() ; ++i) {
  	    if (_names.empty()) {
      		o << "seq"<<i;
	      } else {
      		o << _names[i];
	      }
        for(unsigned j = 0 ; j < this->numsites() ; ++j) {
          o << '\t' << (*this)[i][j];
        }
        if (i < this->size()-1)
          o << '\n';
      }
    }
    else {
    	o << _names[0]; // anc => _names[0]
      for(unsigned i = 0 ; i < this->numsites() ; ++i) {
        o << '\t' << 'N';
      }
    	o << '\n';
      for(unsigned i = 0 ; i < this->size() ; ++i) {
  	    if (_names.empty()) {
      		o << "seq"<<i;
	      } else {
	        o << _names[i+1]; //i => i+1
        }
        for(unsigned j = 0 ; j < this->numsites() ; ++j) {
          o << '\t' << (*this)[i][j];
        }
        if (i < this->size()-1)
          o << '\n';
      }
    }
    return o;
  }
  std::string SimpleSNP::label(unsigned i) const
  { //i is the index in the 'data' matrix
    assert(i<PolyTable::GetData().size());
    return _names[i+1-unsigned(haveOutgroup)];
  }
  void SimpleSNP::set_label(unsigned i, const string & str)
  {
    assert(i<PolyTable::GetData().size());
    _names[i+1-unsigned(haveOutgroup)] = str;
  }
  bool SimpleSNP::outgroup(void) const
  {
    return haveOutgroup;
  }
  void SimpleSNP::RemoveMono(bool skipOutgroup, unsigned outgroup)
  {
    std::vector<double> poss = PolyTable::GetPositions();
    std::vector<std::string> data = PolyTable::GetData();
    std::vector<double> poss_n;
    std::vector<std::string> data_n(data.size());
    for(unsigned j = 0; j < poss.size(); ++ j) {
      stateCounter Counts;
      for(unsigned i = 0; i < data.size(); ++ i) {
        if( (skipOutgroup == false) || (skipOutgroup == true && i != outgroup) )
          Counts( data[i][j] );
      }
      if( Counts.nStates() > 1 ) {
        poss_n.push_back(poss[j]);
        for(unsigned i = 0; i < data.size(); ++ i) 
          data_n[i] += data[i][j];
      }
    }
    assign(&poss_n[0], poss_n.size(), &data_n[0], data_n.size());
  }
  void SimpleSNP::FilterMissing(const int& co_mis, bool skipOutgroup, unsigned outgroup)
  {
    std::vector<double> poss = PolyTable::GetPositions();
    std::vector<std::string> data = PolyTable::GetData();
    std::vector<double> poss_n;
    std::vector<std::string> data_n(data.size());
    for(unsigned j = 0; j < poss.size(); ++ j) {
      stateCounter Counts;
      for(unsigned i = 0; i < data.size(); ++ i) {
        if( (skipOutgroup == false) || (skipOutgroup == true && i != outgroup) )
          Counts( data[i][j] );
      }
      if( (int)Counts.n <= co_mis ) {
        poss_n.push_back(poss[j]);
        for(unsigned i = 0; i < data.size(); ++ i) 
          data_n[i] += data[i][j];
      }
    }
    assign(&poss_n[0], poss_n.size(), &data_n[0], data_n.size());
  }
  void SimpleSNP::subset(IntVec& idxs_ind, IntVec& idxs_pos)
  {
    std::vector<double> poss = PolyTable::GetPositions();
    std::vector<std::string> data = PolyTable::GetData();
    
    if(idxs_ind.size() == 0) push_back(idxs_ind).repeat_fun(data.size(), next<int>(0));
    if(idxs_pos.size() == 0) push_back(idxs_pos).repeat_fun(poss.size(), next<int>(0));
    unsigned n_ind = idxs_ind.size(), n_pos = idxs_pos.size();
    
    vector<double> poss_n;
    for(unsigned j = 0; j < n_pos; ++ j)
      poss_n.push_back( poss[idxs_pos[j]] );
    
    vector<string> labels_n(n_ind+1), data_n(n_ind+haveOutgroup);
    labels_n[0] = _names[0];
    if(haveOutgroup) 
      data_n[0] = data[0];
    for(unsigned i = 0; i < n_ind; ++ i) {
      labels_n[i+1] = _names[idxs_ind[i]+1];
      for(unsigned j = 0; j < n_pos; ++ j)
        data_n[i+haveOutgroup] += data[idxs_ind[i]+haveOutgroup][idxs_pos[j]];
    }

    assign(&poss_n[0], poss_n.size(), &data_n[0], data_n.size());
    _names = labels_n;
  }
  void SimpleSNP::set_outgroup(const string & ind_og)
  {
    std::vector<double> poss = PolyTable::GetPositions();
    std::vector<std::string> data = PolyTable::GetData();
    
    int idx_og = -1;
    for(unsigned i = 0; i < data.size(); i++) {
      if(_names[i+1-haveOutgroup] == ind_og) {
        idx_og = i;
        break;
      }
    }
    if(idx_og < 0) {
      cerr << format("cannot find outgroup: %s\n") % ind_og;
      exit(1);
    }
    
    string data_og = data[idx_og];
    if( haveOutgroup ) {
      _names.erase(_names.begin() + idx_og);
      _names.insert(_names.begin(), ind_og);
      data.erase(data.begin() + idx_og);
      data.insert(data.begin(), data_og);
    } else {
      _names.erase(_names.begin() + idx_og + 1);
      _names[0] = ind_og;
      data.erase(data.begin() + idx_og);
      data.insert(data.begin(), data_og);
      haveOutgroup = true;
    }
    assign(&poss[0], poss.size(), &data[0], data.size());
  }
  void SimpleSNP::clean_outgroup(unsigned & n_N, unsigned & n_fix)
  {
    std::vector<double> poss = PolyTable::GetPositions();
    std::vector<std::string> data = PolyTable::GetData();

    std::vector<double> poss_n;
    std::vector<std::string> data_n(data.size());
    for(unsigned j = 0; j < poss.size(); ++ j) {
      char anc = data[0][j];
      if(anc == 'N') {
        n_N ++;
      } else {
        unsigned n_anc = 0;
        for(unsigned i = 1; i < data.size(); ++ i) {
          if(data[i][j] == anc) n_anc ++;
        }
        if( n_anc == 0 ) {
          n_fix ++;
        } else {
          poss_n.push_back(poss[j]);
          for(unsigned i = 0; i < data.size(); ++ i) 
            data_n[i] += data[i][j];
        }
      }
    }
    assign(&poss_n[0], poss_n.size(), &data_n[0], data_n.size());
  }
  void SimpleSNP::write(const string& fo, const string& out_format, const string& chr)
  {
    std::vector<double> poss = PolyTable::GetPositions();
    std::vector<std::string> data = PolyTable::GetData();
    ofstream fho(fo.c_str());
    unsigned n_ind = data.size(), n_pos = poss.size();

    if(out_format == "ssp") 
    {
      print(fho);
    }
    else if(out_format == "phylip")
    {
      fho << format("%d  %d\n") % n_ind % n_pos;
      for (uint32_t r = 0; r < ceil(double(n_pos) / 250); ++r) {
        for (uint32_t i = 0; i < n_ind; i++) {
          string data_row_string = "";
          for (size_t j = 250*r; j < min(n_pos, 250*r+250); j++)
            data_row_string += data[i][j];
          if(r == 0) {
            fho << format("%-10s%s\n") % label(i) % data_row_string;
          } else {
            fho << format("%s\n") % data_row_string;
          }
        }
        fho << "\n";
      }
    }
    else if(out_format == "table")
    {
      for (PolyTable::const_site_iterator it=sbegin(); it!=send(); ++it) {
        double pos = it->first;
        string str = it->second;
        string str_new = str.substr(str.length()-1, 1);
        for(uint32_t k = 0; k < str.length() - 1; k ++)
          str_new += "\t" + str.substr(k, 1);
        fho << format("%.0f\t%s\n") % pos % str_new;
      }
    }
    else if(out_format == "structure")
    {
      string row1;
      for(unsigned j = 0; j < n_pos; j++) 
        row1 += boost::lexical_cast<string>(poss[j]) + "\t";
      fho << format("%s\n") % row1;
      for(unsigned i = 0; i < n_ind; i ++) {
        string row = str(format("%s\t%s\t") % label(i) % (i+1)) ;
        for(unsigned j = 0; j < n_pos; j ++) {
          switch (char(std::toupper(data[i][j]))) {
            case 'A':
              row += "1\t";
              break;
            case 'T':
              row += "2\t";
              break;
            case 'C':
              row += "3\t";
              break;
            case 'G':
              row += "4\t";
              break;
            case 'N':
              row += "9\t";
              break;
            default:
              row += "9\t";
          }
        }
        fho << format("%s\n%s\n") % row % row;
      }
    }
    else if(out_format == "ped")
    {
      string fo2 = fo + ".map";
      ofstream fho2(fo2.c_str());
      for(uint32_t j = 0; j < n_pos; j++)
        fho2 << format("%s\t%d\t0\t%d\n") % chr % (j+1) % poss[j];
      for(unsigned i = 0; i < n_ind; i++) {
        string id = label(i);
        string row = id + "\t" + boost::lexical_cast<string>(i+1) + "\t0\t0\t0\t0\t";
        for(unsigned j = 0; j < n_pos; j ++) {
          string nt = boost::lexical_cast<string>(data[i][j]);
          row += nt + " " + nt + "\t";
        }
        fho << format("%s\n") % row;
      }
    }
    else
    {
      cerr << format("unsupported format: %s\n") % out_format;
      exit(1);
    }
    fho.close();
  }
  void SimpleSNP::write_expand(const string& dirO, const string& seqid)
  {
    std::vector<double> poss = PolyTable::GetPositions();
    std::vector<std::string> data = PolyTable::GetData();
    uint32_t n_ind = data.size(), n_pos = poss.size();

    for(unsigned i = 0; i < n_ind; i++) {
      string id = label(i+1);
      string fo = dirO + "/" + id + "/" + seqid;
      
      ofstream fho(fo.c_str());
      if (!fho.is_open()) {
        cout << format("cannot write output: %s\n") % fo;
        exit(1);
      }
      for (unsigned j = 0; j < n_pos; j ++) {
        fho << format("%d\t%s\n") % int(poss[j]) % data[i][j];
      }
      fho.close();
    }
  }
}



