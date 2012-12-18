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
#include "ssp.h"
using boost::format;
using namespace std;
using namespace boost::assign;
using namespace Sequence;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string fi, fo, out_format, chr;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (simpleSNP) file")
    ("out,o", po::value<string>(&fo), "output file")
    ("out_format,f", po::value<string>(&out_format)->default_value("phylip"), "supported output format: plain, ssp, phylip, table, structure, ped")
    ("chr,c", po::value<string>(&chr)->default_value("chr0"), "chromosome (required by PED format)")
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
  if(out_format == "plain") {
    vector<string> ss;
    boost::split(ss, fi, boost::is_any_of("/"));
    boost::split(ss, ss[ss.size()-1], boost::is_any_of("."));
    ssp.write_expand(fo, ss[0]);
  } else {
    ssp.write(fo, out_format, chr);
  }
 
  return EXIT_SUCCESS;
}



