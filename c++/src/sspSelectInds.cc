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

void remove_mono(SimpleSNP& ssp) {
  uint32_t n_ind = ssp.size(), n_pos = ssp.numsites();
  vector<string> data = ssp.GetData();
  vector<double> poss = ssp.GetPositions();
  vector<string> data_n(n_ind);
  vector<double> poss_n;
  for(uint32_t j = 0; j < n_pos; j++) {
    set<string> nts;
    bool tag_mono = true;
    for(uint32_t i = 0; i < n_ind; i++) {
      string nt = boost::lexical_cast<string>(data[i][j]);
      if(nt != "N") {
        nts.insert( nt );
        if(nts.size() > 1) {
          tag_mono = false;
          break;
        }
      }
    }
    if( tag_mono == false ) {
      poss_n.push_back(poss[j]);
      for(uint32_t i = 0; i < n_ind; i++)
        data_n[i] += data[i][j];
    }
  }
  if( !ssp.assign(&poss_n[0], poss_n.size(), &data_n[0], data_n.size()) ) {
    cout << "error assigning data to ssp\n";
    exit(1);
  }
}

int main(int argc, char *argv[]) {
  string f_in, f_out, f_acc, opt_ind;
  int opt_conf;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("f_in,i", po::value<string>(&f_in), "input (simpleSNP) file")
    ("f_out,o", po::value<string>(&f_out), "output (simpleSNP) file")
    ("f_acc,a", po::value<string>(&f_acc)->default_value("/project/youngn/zhoup/Scripts/conf/acc_ids.tbl"), "acc option file")
    ("opt_ind,t", po::value<string>(&opt_ind), "id option")
    ("opt_conf,c", po::value<int>(&opt_conf)->default_value(1), "id_all option")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("f_in") || !vm.count("f_out") || !vm.count("opt_ind")) {
    cout << cmdOpts << endl;
    return 1;
  }

  SimpleSNP ssp;
  ifstream fh01( f_in.c_str() );
  ssp.read( fh01 );
  fh01.close();

  uint32_t n_ind = ssp.size(), n_pos = ssp.numsites();
  cout << format("before ind_select:\n  %d inds, %d sites\n") % n_ind % n_pos;
 
  IntVec idxs_pos, idxs_ind;
  StrVec inds;
  idxs_ind = get_acc_idx(f_acc, opt_ind, inds, opt_conf);
  
  ssp.subset(idxs_ind, idxs_pos);
  ssp.RemoveMono();
  vector<string> data = ssp.GetData();
  vector<double> poss = ssp.GetPositions();
  n_ind = ssp.size(), n_pos = ssp.numsites();

  cout << format("after ind_select:\n  %d inds, %d sites\n") % n_ind % n_pos;

  ofstream fho(f_out.c_str());
  ssp.print(fho);
  /*
  fho << format("%d\t%d\n") % n_ind % n_pos;
  string row1 = "", row2 = "";
  for(uint32_t i = 0; i < n_pos; i ++) {
    row1 += boost::lexical_cast<string>(poss[i]) + "\t";
    row2 += "N\t";
  }
  fho << format("%s\n%s\n") % row1 % row2;
  for(uint32_t i = 0; i < n_ind; i++) {
    string row = inds[i] + "\t";
    for (uint32_t j = 0; j < n_pos; j ++)
      row += boost::lexical_cast<string>(data[i][j]) + "\t";
    fho << format("%s\n") % row;
  }*/
  fho.close();

  return EXIT_SUCCESS;
}



