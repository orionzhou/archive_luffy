#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <map>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "bam_utils.h" 
using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
using namespace BamTools;


int main(int argc, char *argv[]) {
  string fi, fo, regionStr, tag;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("f_in,i", po::value<string>(&fi), "input (BAM) file")
    ("f_out,o", po::value<string>(&fo), "output file")
    ("region,r", po::value<string>(&regionStr), "region")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("f_in") || !vm.count("f_out") || !vm.count("region")) {
    cout << cmdOpts << endl;
    return 1;
  }
  
  Utilities utils;
  BamReader reader;
  BamRegion region;
  BamWriter writer;
  
  reader.Open(fi);
  utils.ParseRegionString(regionStr, reader, region); 
  reader.SetRegion(region);
  cout << format("region set to: %s\n") % regionStr;
  ofstream fho( fo.c_str() );

  map<string, int> rm;
  map<string, int>::iterator rm_it;
  pair< map<string, int>::iterator, bool > rm_p;

  BamAlignment al;
  vector<CigarOp> cigars;
  int refId = region.LeftRefID;
  while( reader.GetNextAlignment(al) ) {
    if( al.IsDuplicate() || al.MapQuality == 0) continue;
    cigars =  al.CigarData;
    int type = 0;
// type = 1(unmapped) 2(mapped on a different chr) 3(soft-clipped)
    if( al.IsMapped() && !al.IsMateMapped() ) {
      type = 1;
    } else if( al.IsMapped() && al.IsMateMapped() ) {
      if ( al.MateRefID != refId ) {
        type = 2;
      } else if( (cigars[0].Type == 'S' && cigars[0].Length >= 10) ||
        (cigars[cigars.size()-1].Type == 'S' && cigars[cigars.size()-1].Length >= 10) ) {
        type = 3;
      }
    }

    if(type > 0) {
      rm_p = rm.insert( pair<string, int> (al.Name, type) );
      if(rm_p.second == false) {
        cout << format("reading %s twice\n") % al.Name;
        exit(0);
      }
      int orphanIsFirst = al.IsFirstMate() ? 0 : 1;
      if(type == 3) orphanIsFirst = al.IsFirstMate() ? 1 : 0; 
      fho << format("%s\t%d\t%d\n") % al.Name % type % orphanIsFirst;
    }
  }
  reader.Close();

  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
