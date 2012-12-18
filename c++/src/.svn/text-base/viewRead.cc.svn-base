#include <iostream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "bam_utils.h"
using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
using namespace BamTools;

int main(int argc, char *argv[]) {
  string fi, fo, dirW, regionStr, name;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (BAM) file")
    ("region,r", po::value<string>(&regionStr), "region")
    ("name,n", po::value<string>(&name)->default_value(""), "read name")
    ("out,o", po::value<string>(&fo)->default_value("tmp.bam"), "output file")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("name")) {
    cout << cmdOpts << endl;
    return 1;
  }

  fs::path fi1( fi );
  fs::path fi2( fi + ".bai" );
  BamReader reader;
  BamWriter writer;
  Utilities util;
  if( exists(fi1) ) {
    reader.Open(fi1.string());
    if( reader.HasIndex() ) {
      cout << format("%s opened with index ...\n") % fi;
      if(regionStr != "") {
        BamRegion region;
        if( util.ParseRegionString(regionStr, reader, region) ) {
          reader.SetRegion(region);
          cout << format("region set to: %s\n") % regionStr;
        } else {
          cerr << format("invalid region string: %s\n") % regionStr;
          exit(1);
        }
      } else {
        cout << "no region specified: starting from scratch\n";
      }
    } else {
      cout << format("%s opened without index ...\n") % fi;
    }
  } else {
    cout << format("cannot read from %s ...\n") % fi;
    exit(1);
  }
  writer.Open( fo, reader.GetHeader(), reader.GetReferenceData() );

  BamAlignment al;
  bool tag_s = false, tag_e = false;
  uint32_t n_more = 0, cnt = 0, nout = 10;
  while( reader.GetNextAlignment(al) ) {
    string name2(al.Name);
    if(name == name2) tag_s = true;
    if(tag_s == true) {
      if(name != name2) tag_e = true;
      util.PrintRead(al);
      writer.SaveAlignment(al);
    }
    cnt ++;
    if(tag_e == true && ++n_more >= nout) break;
  }
  reader.Close();
  writer.Close();
  cout << format("%10.0f alignments read\n") % cnt;
  cout << format("%10.0f alignments written to %s\n") % nout % fo;
  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}


 
