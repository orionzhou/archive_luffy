#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include "ssp.h" 
using boost::format;
using namespace std;
using namespace boost::assign;
using namespace Sequence;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    string fi, fo, fa, opt;
    po::options_description cmdOpts("Allowed options");
    cmdOpts.add_options()
        ("help,h", "produce help message")
        ("fi,i", po::value<string>(&fi), "input (simpleSNP) file")
        ("fo,o", po::value<string>(&fo), "output (simpleSNP) file")
        ("fa,a", po::value<string>(&fa)->default_value("/home/youngn/zhoup/git/conf/acc_ids.tbl"), "acc option file")
        ("opt,p", po::value<string>(&opt_ind), "id option")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
    po::notify(vm);
    if(vm.count("help") || !vm.count("fin") || !vm.count("fo") || !vm.count("opt")) {
        cout << cmdOpts << endl;
        return 1;
    }

    istream* in = (!vm.count("in") || fi == "" || fi == "-" || fi == "stdin") ?
        &cin : new ifstream( fi.c_str() );
    ostream* out = (!vm.count("out") || fo == "" || fo == "-" || fo == "stdout") ?
        &cout : new ofstream( fo.c_str() );
   
    SSP ssp;
    ssp.read( *in );
    cerr << format("before ind_select:\n  %d inds, %d sites\n") % ssp.nind % ssp.npos;
  
    IntVec idxs_pos, idxs_ind;
    StrVec inds_a = get_all_labels(ssp);
    cerr << format("  labels: %s\n") % boost::join(inds_a, " ");
    StrVec inds = get_acc_ids(fa, opt);
    idxs_ind = get_idx(inds, inds_a);
    
    ssp.subset(idxs_ind, idxs_pos);
    ssp.RemoveMono();
    cerr << format("after ind_select:\n  %d inds, %d sites\n") % ssp.nind % ssp.npos;

    ssp.write(*out, "ssp");

    return EXIT_SUCCESS;
}



