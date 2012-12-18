#include <iostream>
#include <stdio.h>
#include <vector>
#include <map>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/PolyTableSlice.hpp>
#include "ssp.h"

using namespace std;
using boost::format;
using namespace Sequence;
using namespace Sequence::Alignment;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string fi, fo;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input file")
    ("out,o", po::value<string>(&fo), "output file")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);

  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  fs::path f01( fi );
  if ( !exists( f01 ) ) {
    cout << "\nNot found: " << f01.string() << endl;
    return 1;
  } else {
    cout << format("reading from %s\n") % fi;
  }

  SimpleSNP ssp;
  ifstream fh01( f01.string().c_str() );
  ssp.read( fh01 );
  fh01.close();
  cout << format("\t%d positions, %d sequences\n") % ssp.numsites() % ssp.size();
  PolySites psi(ssp.sbegin(), ssp.send());
//  PolySNP psp(&psi);
//  cout << "ThetaPi: " << psp.ThetaPi() << endl;
//  cout << "ThetaW: " << psp.ThetaW() << endl;
//  cout << "Tajima's D: " << psp.TajimasD() << endl;
  double chrLen = *(psi.pbegin()+psi.numsites()-1) + 1;
  unsigned winSize = 100000; 
  unsigned winStep = 100000; 
  PolyTableSlice<PolySites> pts(psi.sbegin(), psi.send(), winSize, winStep, chrLen);

  ofstream fho(fo.c_str());
  if (!fho.is_open()) {
    cout << format("cannot write output: %s\n") % fo;
    return 1;
  }

  for (unsigned i = 0; i < pts.size(); i++) {
    PolySites psi1(pts[i]);
    //PolySNP psp1(&psi);
    //double thetaW = winCalc.ThetaW() / winSize;
    if(psi1.numsites() > 0) {
      fho << format("%.0f\t%.0f\n") % *(psi1.pbegin()) % *(psi1.pbegin()+psi1.numsites()-1);
    } else {
      fho << format("\t\n");
    }
  }
  fho.close();
  cout << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
	return EXIT_SUCCESS;
}


