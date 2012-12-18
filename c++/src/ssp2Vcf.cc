#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>
#include <Sequence/stateCounter.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/PolyTableSlice.hpp>
#include "ssp.h"
using boost::format;
using namespace std;
using namespace boost::assign;
using namespace Sequence;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string fi, fo, ind_ref, chr;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (simpleSNP) file")
    ("out,o", po::value<string>(&fo), "output (VCF) file")
    ("ref,r", po::value<string>(&ind_ref)->default_value(""), "reference ind")
    ("chr,c", po::value<string>(&chr)->default_value("chrU"), "chr")
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

  uint32_t n_ind = ssp.size(), n_pos = ssp.numsites();
  cout << format("%s:\n\t%d inds, %d sites\n") % fi % n_ind % n_pos;
 
  int idx_ref = -1;
  if(ind_ref == "") ind_ref = ssp.label(1);
  for(unsigned i = 0; i < n_ind; i++) {
    string ind = ssp.label(i+1);
    if(ind == ind_ref) {
      idx_ref = i;
      break;
    }
  }
  if(idx_ref < 0) {
    cout << format("cannot find reference: %s\n") % ind_ref;
    exit(1);
  }

  ofstream fho( fo.c_str() );
  vector<string> data = ssp.GetData();
  vector<double> poss = ssp.GetPositions();
  for(unsigned j = 0; j < poss.size(); ++ j) {
    char ref = data[idx_ref][j];
    if(ref == 'N') {
      cout << format("position %d skipped : ref = N\n") % poss[j];
      continue;
    }
    stateCounter Counts;
    for(unsigned i = 0; i < data.size(); ++ i) Counts( data[i][j] );
    if(ref != 'A' && Counts.a > 0) fho << format("%s\t%.0f\t%s\t%s\t%s\n") % chr % poss[j] % "" % ref % "A";
    if(ref != 'T' && Counts.t > 0) fho << format("%s\t%.0f\t%s\t%s\t%s\n") % chr % poss[j] % "" % ref % "T";
    if(ref != 'C' && Counts.c > 0) fho << format("%s\t%.0f\t%s\t%s\t%s\n") % chr % poss[j] % "" % ref % "C";
    if(ref != 'G' && Counts.g > 0) fho << format("%s\t%.0f\t%s\t%s\t%s\n") % chr % poss[j] % "" % ref % "G";
  }

  return EXIT_SUCCESS;
}



