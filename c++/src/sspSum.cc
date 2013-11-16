#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include "ssp.h" 
using namespace std;
using boost::format;
using namespace Sequence;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
    string fi, fo;
    po::options_description cmdOpts("Allowed options");
    cmdOpts.add_options()
        ("help,h", "produce help message")
        ("in,i", po::value<string>(&fi)->implicit_value(""), "input (SimpleSNP format)")
        ("out,o", po::value<string>(&fo)->implicit_value(""), "output")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
    po::notify(vm);

    if(vm.count("help")) {
        cout << cmdOpts << endl;
        return 1;
    }
    
    istream* in = (!vm.count("in") || fi == "" || fi == "-" || fi == "stdin") ?
        &cin : new ifstream( fi.c_str() );
    ostream* out = (!vm.count("out") || fo == "" || fo == "-" || fo == "stdout") ?
        &cout : new ofstream( fo.c_str() );
   
    SSP ssp;
    ssp.read( *in );
    
    bool hasOG = ssp.outgroup();
    *out << "Outgroup: " << (hasOG ? "Yes" : "No") << endl;
    *out << ssp.npos << " sites" << endl;
    *out << ssp.nind << " samples: " << boost::algorithm::join(ssp.labels, ",") << endl;
    return EXIT_SUCCESS;
}
