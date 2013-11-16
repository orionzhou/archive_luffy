#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <Sequence/stateCounter.hpp>
#include "ssp.h" 
using namespace std;
using boost::format;
using namespace Sequence;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    string fi, fo;
    po::options_description cmdOpts("Allowed options");
    cmdOpts.add_options()
        ("help,h", "produce help message")
        ("in,i", po::value<string>(&fi)->implicit_value(""), "input (SimpleSNP format)")
        ("out,o", po::value<string>(&fo)->implicit_value(""), "output")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
    po::notify(vm);

    if(vm.count("help")) {
        cout << cmdOpts << endl;
        return 1;
    }
    
    istream* in = (!vm.count("in") || fi == "" || fi == "-" || fi == "stdin") ?
        &cin : new ifstream( fi.c_str() );
    ostream* out = (!vm.count("out") || fo == "" || fo == "-" || fo == "stdout") ?
        &cout : new ofstream( fo.c_str() );
   
    SSP ssp;
    ssp.read( *in );
    
    vector<string> data = ssp.GetData();
    vector<double> poss = ssp.GetPositions();
    bool hasOG = ssp.outgroup();
    unsigned nind = ssp.size();
    unsigned npos = poss.size();
    
    *out << "pos\tn_states\tanc\tder\tn_N\tn_anc\tn_der" << endl;
    for(unsigned j = 0; j < npos; ++j) {
        stateCounter Counts;
        for(unsigned i = 0; i < nind; ++i) Counts( data[i][j] );
        
        int n_states = Counts.nStates();
        int n_N = Counts.n;
        if(n_states <= 1) continue;

        map<char, int> am;
        map<char, int>::iterator am_it;
        if(Counts.a > 0) am['A'] = Counts.a;
        if(Counts.t > 0) am['T'] = Counts.t;
        if(Counts.c > 0) am['C'] = Counts.c;
        if(Counts.g > 0) am['G'] = Counts.g;
        
        char anc, der, allele;
        int n_anc=0, n_der=0, count;
        if(hasOG) {
            anc = char(std::toupper(data[0][j]));
            if(anc == 'N') continue;
            for(am_it=am.begin(); am_it != am.end(); am_it++) {
                allele = am_it->first;
                count = am_it->second;
                if(allele == anc) {
                    n_anc = count - 1;
                } else {
                    if(n_der < count) {
                        der = allele;
                        n_der = count;
                    }
                }
            }
        } else {
            for(am_it=am.begin(); am_it != am.end(); am_it++) {
                allele = am_it->first;
                count = am_it->second;
                if(n_anc < count) {
                    anc = allele;
                    n_anc = count;
                }
            }
            for(am_it=am.begin(); am_it != am.end(); am_it++) {
                allele = am_it->first;
                count = am_it->second;
                if(allele == anc) {
                } else if(n_der < count) {
                    der = allele;
                    n_der = count;
                }
            }
        }
        *out << format("%d\t%d\t%s\t%s\t%d\t%d\t%d") % poss[j] % n_states % anc % der % n_N % n_anc % n_der << endl;
    }
    return EXIT_SUCCESS;
}
