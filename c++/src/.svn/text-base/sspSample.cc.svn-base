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
  string fi, fo, out_format;
  uint32_t n_sam;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (simpleSNP) file")
    ("out,o", po::value<string>(&fo), "output file")
    ("out_format,f", po::value<string>(&out_format)->default_value("ssp"), "output format")
    ("sample,n", po::value<uint32_t>(&n_sam)->default_value(0), "# positions to sample")
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
  ssp.read(fh01);
  fh01.close();

  unsigned n_ind = ssp.size(), n_pos = ssp.numsites();
  cout << format("before sampling:  %3d inds, %6d positions\n") % n_ind % n_pos;
  n_sam = (n_sam == 0) ? n_pos : n_sam;

  IntVec idxs_ind;
  IntVec idxs_pos = sample_serial(n_pos, n_sam);
  ssp.subset(idxs_ind, idxs_pos);
  
  n_pos = ssp.numsites();
  cout << format("after sampling:  %3d inds, %6d positios\n") % n_ind % n_pos;

  ssp.write(fo, out_format);

  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}



