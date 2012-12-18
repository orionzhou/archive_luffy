#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>
#include <Sequence/stateCounter.hpp>
#include "ssp.h"
using namespace std;
using boost::format;
using namespace Sequence;
using namespace Sequence::Alignment;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string fi, fo;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (simpleSNP) file")
    ("out,o", po::value<string>(&fo), "output file")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);

  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  SimpleSNP ssp;
  ifstream fh01( fi.c_str() );
  ssp.read( fh01 );
  fh01.close();

  ofstream fho(fo.c_str());

  vector<string> data = ssp.GetData();
  vector<double> poss = ssp.GetPositions();
  bool hasOG = ssp.outgroup();

  fho << "pos\tn_states\tanc\tder\tn_N\tn_anc\tn_der\n";
  for(unsigned j = 0; j < poss.size(); ++j) {
    stateCounter Counts;
    for(unsigned i = 0; i < data.size(); ++i) Counts( data[i][j] );
    
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
    fho << format("%d\t%d\t%s\t%s\t%d\t%d\t%d\n") % poss[j] % n_states % anc % der % n_N % n_anc % n_der;
  }
  fho.close();
	return EXIT_SUCCESS;
}
