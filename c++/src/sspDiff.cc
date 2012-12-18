#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>
#include "ssp.h"
using namespace std;
using boost::format;
using namespace boost::assign;
using namespace Sequence;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string fi, fo, out_format, ind1, ind2;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (simpleSNP) file")
//    ("out,o", po::value<string>(&fo), "output file")
//    ("out_format,f", po::value<string>(&out_format)->default_value("ssp"), "output format")
    ("ind1,1", po::value<string>(&ind1), "individual 1")
    ("ind2,2", po::value<string>(&ind2), "individual 2")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("ind1") || !vm.count("ind2")) {
    cout << cmdOpts << endl;
    return 1;
  }

  SimpleSNP ssp;
  ifstream fh01( fi.c_str() );
  ssp.read(fh01);
  fh01.close();

  uint32_t n_ind = ssp.size(), n_pos = ssp.numsites();

  int idx1 = -1, idx2 = -1;
  for(int i = 0; i < (int)ssp.size(); i++) {
    string ind = ssp.label(i+1);
    if(ind == ind1) idx1 = i;
    if(ind == ind2) idx2 = i;
  }
  if(idx1 < 0) {
    cerr << format("no individual called %s\n") % ind1;
    exit(1);
  }
  if(idx2 < 0) {
    cerr << format("no individual called %s\n") % ind2;
    exit(1);
  }

  int n_same=0, n_diff=0, n_1n=0, n_2n=0, n_nn=0;
  vector<string> data = ssp.GetData();
  for(uint32_t j = 0; j < n_pos; ++j) {
    char s1 = char(std::toupper(data[idx1][j]));
    char s2 = char(std::toupper(data[idx2][j]));
    if(s1 == 'N' && s2 == 'N') {
      n_nn ++;
    } else if(s1 == 'N' && s2 != 'N') {
      n_1n ++;
    } else if(s1 != 'N' && s2 == 'N') {
      n_2n ++;
    } else if(s1 == s2) {
      n_same ++;
    } else {
      n_diff ++;
    }
  }
  
  cout << format("A total of %3d inds, %6d positions\n") % n_ind % n_pos;
  cout << format("%6d positions: %s  = %s\n") % n_same % ind1 % ind2;
  cout << format("%6d positions: %s != %s\n") % n_diff % ind1 % ind2;
  cout << format("%6d positions: %s  = N && %s != N\n") % n_1n % ind1 % ind2;
  cout << format("%6d positions: %s != N && %s  = N\n") % n_2n % ind1 % ind2;
  cout << format("%6d positions: %s  = N && %s  = N\n") % n_nn % ind1 % ind2;

  return EXIT_SUCCESS;
}



