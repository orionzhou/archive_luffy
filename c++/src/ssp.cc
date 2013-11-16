#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <cctype>
#include <cassert>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/stateCounter.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/lexical_cast.hpp>
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
IntVec get_idx(StrVec& svi, StrVec& sva) {
        cout << format("  selecting ids: %s\n") % boost::join(svi, " ");
        map<string, int> tmp;
        for(uint32_t i = 0; i < sva.size(); i ++) tmp[sva[i]] = i;

        IntVec idxs;
        for(uint32_t i = 0; i < svi.size(); i ++) {
                string id = svi[i];
                map<string, int>::iterator it = tmp.find(id);
                if(it == tmp.end()) {
                        cout << format("  %s not found\n") % id;
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
    std::istream & SSP::read (std::istream& s) throw (Sequence::badFormat, std::exception) {
        SimpleSNP::read(s);
        nind = SimpleSNP::size();
        npos = SimpleSNP::numsites();
        for(unsigned i = 0; i < SimpleSNP::size(); i ++)
            labels.push_back( SimpleSNP::label(i+1) );
        return s;
    }
    std::ostream& SSP::write (std::ostream& o, const string& fmt, const string& chr) {
        bool haveOG = SimpleSNP::outgroup();
        std::vector<double> poss = SimpleSNP::GetPositions();
        std::vector<std::string> data = SimpleSNP::GetData();
        if(fmt == "ssp") {
            o << nind-haveOG << '\t' << npos << endl;
            
            for(unsigned j = 0; j < npos-1; j++)
                o << SimpleSNP::position(j) << '\t';
            o << SimpleSNP::position(npos-1) << endl;
            
            if(!haveOG) {
                o << "anc" << '\t';
                for(unsigned j = 0; j < npos-1; j++)
                    o << "N" << '\t';
                o << "N" << endl;
            }
            for(unsigned i = 0; i < nind; i ++) {
                o << labels[i] << '\t';
                for(unsigned j = 0; j < npos-1; j ++)
                    o << data[i][j] << '\t';
                o << data[i][npos-1] << endl;
            }
        } else if(fmt == "phylip") {
            o << format("%d  %d\n") % nind % npos;
            for (uint32_t r = 0; r < ceil(double(npos) / 250); ++r) {
                for (uint32_t i = 0; i < nind; i++) {
                    string data_row_string = "";
                    for (size_t j = 250*r; j < min(npos, 250*r+250); j++)
                        data_row_string += data[i][j];
                    if(r == 0) {
                        o << format("%-10s%s\n") % labels[i] % data_row_string;
                    } else {
                        o << format("%s\n") % data_row_string;
                    }
                }
                o << endl;
            }
        } else if(fmt == "table") {
            for (PolyTable::const_site_iterator it=sbegin(); it!=send(); ++it) {
                double pos = it->first;
                string str = it->second;
                string str_new = str.substr(str.length()-1, 1);
                for(uint32_t k = 0; k < str.length() - 1; k ++)
                    str_new += "\t" + str.substr(k, 1);
                o << format("%.0f\t%s") % pos % str_new << endl;
            }
        } else if(fmt == "structure") {
            string row1;
            for(unsigned j = 0; j < npos; j++) 
                row1 += boost::lexical_cast<string>(poss[j]) + "\t";
            o << format("%s\n") % row1;
            for(unsigned i = 0; i < nind; i ++) {
                string row = str(format("%s\t%s\t") % labels[i] % (i+1)) ;
                for(unsigned j = 0; j < npos; j ++) {
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
                o << format("%s\n%s") % row % row << endl;
            }
        } else if(fmt == "ped") {
            for(uint32_t j = 0; j < npos; j++)
                o << format("%s\t%d\t0\t%d\n") % chr % (j+1) % poss[j];
            for(unsigned i = 0; i < nind; i++) {
                string id = labels[i];
                string row = id + "\t" + boost::lexical_cast<string>(i+1) + "\t0\t0\t0\t0\t";
                for(unsigned j = 0; j < npos; j ++) {
                    string nt = boost::lexical_cast<string>(data[i][j]);
                    row += nt + " " + nt + "\t";
                }
                o << row << endl;
            }
        } else {
            cerr << format("unsupported format: %s\n") % fmt;
            exit(1);
        }
        return o;
    }
    void SSP::set_label(unsigned i, const string& label) {
        labels[i-1] = label;
    }
    unsigned SSP::get_label_idx(const string& label) {
        unsigned idx = -1;
        for(unsigned i = 0; i < nind; i++) {
            if(labels[i] == label) {
                idx = i;
                break;
            }
        }
        if(idx < 0) {
            cerr << "label not find: " << label << endl;
            exit(1);
        }
        return idx;
    }

    void SSP::ApplyFreqFilter(unsigned mincount,bool haveOutgroup, unsigned outgroup) {
        SimpleSNP::ApplyFreqFilter(mincount, haveOutgroup, outgroup);
        npos = SimpleSNP::numsites();
    }
    void SSP::RemoveMultiHits(bool skipOutgroup, unsigned outgroup) {
        SimpleSNP::RemoveMultiHits(skipOutgroup, outgroup);
        npos = SimpleSNP::numsites();
    }
    void SSP::RemoveMissing(bool skipOutgroup, unsigned outgroup) {
        SimpleSNP::RemoveMissing(skipOutgroup, outgroup);
        npos = SimpleSNP::numsites();
    }
    void SSP::RemoveAmbiguous(bool skipOutgroup, unsigned outgroup) {
        SimpleSNP::RemoveAmbiguous(skipOutgroup, outgroup);
        npos = SimpleSNP::numsites();
    }

    void SSP::FilterMissing(const int& co_mis, bool skipOutgroup, unsigned outgroup) {
        std::vector<double> poss = SimpleSNP::GetPositions();
        std::vector<std::string> data = SimpleSNP::GetData();
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
        SimpleSNP::assign(&poss_n[0], poss_n.size(), &data_n[0], data_n.size());
        npos = SimpleSNP::numsites();
    }
    void SSP::RemoveMono(bool skipOutgroup, unsigned outgroup) {
        std::vector<double> poss = SimpleSNP::GetPositions();
        std::vector<std::string> data = SimpleSNP::GetData();
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
        SimpleSNP::assign(&poss_n[0], poss_n.size(), &data_n[0], data_n.size());
        npos = SimpleSNP::numsites();
    }

    void SSP::subset(IntVec& idxs_ind, IntVec& idxs_pos) {
        bool haveOG = SimpleSNP::outgroup();
        std::vector<double> poss = PolyTable::GetPositions();
        std::vector<std::string> data = PolyTable::GetData();
        
        if(idxs_ind.size() == 0) push_back(idxs_ind).repeat_fun(data.size(), next<int>(0));
        if(idxs_pos.size() == 0) push_back(idxs_pos).repeat_fun(poss.size(), next<int>(0));
        unsigned n_ind = idxs_ind.size(), n_pos = idxs_pos.size();
        
        vector<double> poss_n;
        for(unsigned j = 0; j < n_pos; ++ j)
            poss_n.push_back( poss[idxs_pos[j]] );
        
        vector<string> labels_n(n_ind+haveOG), data_n(n_ind+haveOG);
        if(haveOG) {
            labels_n[0] = labels[0];
            data_n[0] = data[0];
        }
        for(unsigned i = 0; i < n_ind; ++ i) {
            labels_n[i+haveOG] = labels[idxs_ind[i]+haveOG];
            for(unsigned j = 0; j < n_pos; ++ j)
                data_n[i+haveOG] += data[idxs_ind[i]+haveOG][idxs_pos[j]];
        }

        SimpleSNP::assign(&poss_n[0], poss_n.size(), &data_n[0], data_n.size());
        labels = labels_n;
        nind = SimpleSNP::size();
        npos = SimpleSNP::numsites();
    }

/*    void SimpleSNP::write_expand(const string& dir, const string& seqid) {
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
*/
}



