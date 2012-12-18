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
  string fi, fo, ind_og;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (simpleSNP) file")
    ("out,o", po::value<string>(&fo), "output (simpleSNP) file")
    ("outgroup,g", po::value<string>(&ind_og), "outgroup")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("out") || !vm.count("outgroup")) {
    cout << cmdOpts << endl;
    return 1;
  }

  SimpleSNP ssp;
  ifstream fh01( fi.c_str() );
  ssp.read(fh01);
  fh01.close();

  unsigned n_ind = ssp.size(), n_pos = ssp.numsites();
  cout << format("input file: %3d inds, %6d positions, %1d outgroup\n") % n_ind % n_pos % unsigned(ssp.outgroup());

  cout << format("  setting %s to outgroup\n") % ind_og;
  ssp.set_outgroup( ind_og );
  unsigned n_N = 0, n_fix = 0;
  ssp.clean_outgroup(n_N, n_fix);
  cout << format("    %d positions removed: %s = 'N'\n") % n_N % ind_og;
  cout << format("    %d positions removed: fixed in ingroup\n") % n_fix;

  ofstream fho(fo.c_str());
  ssp.print(fho);
  fho.close();

  return EXIT_SUCCESS;
}



