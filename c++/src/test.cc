#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <map>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <api/BamReader.h>
#include <api/BamWriter.h>
using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
using namespace BamTools;

int main( int argc, char* argv[] ) {
  string fi, fo;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi)->default_value("/tmp"), "input file")
    ("out,o", po::value<string>(&fo), "output file")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    return 1;
  }
  
  BamReader reader;
  if(!reader.Open(fi) || !reader.OpenIndex(fi+".bai") || !reader.HasIndex()) {
    cout << "error opening " << fi << endl;
    exit(1);
  }
  
  int refId = reader.GetReferenceID("chr5");
  BamRegion region(refId, 73042, refId, 75094);
  reader.SetRegion(region);
  BamAlignment al;
  reader.GetNextAlignment(al); 
  cout << al.Name << "\t" << al.Position << endl;
 
  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
