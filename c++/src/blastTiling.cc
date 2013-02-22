#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <map>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "common.h"
using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

struct BlastRecord {
    string qId;
    uint32_t qBeg;
    uint32_t qEnd;
    string strand;
    string hId;
    uint32_t hBeg;
    uint32_t hEnd;
    float pct;
    double e;
    float score;
    uint32_t len;
    uint32_t mismatch;
    uint32_t gap;
};
typedef vector<BlastRecord> BlastRecords;

BlastRecord make_blast_record(const vector<string>& ss) {   
    BlastRecord br;
    br.qId = ss[0];
    br.hId = ss[1];
    br.pct = boost::lexical_cast<float>(ss[2]);
    br.len = boost::lexical_cast<uint32_t>(ss[3]);
    br.mismatch = boost::lexical_cast<uint32_t>(ss[4]);
    br.gap = boost::lexical_cast<uint32_t>(ss[5]);
    br.qBeg = boost::lexical_cast<uint32_t>(ss[6]);
    br.qEnd = boost::lexical_cast<uint32_t>(ss[7]);
    br.hBeg = boost::lexical_cast<uint32_t>(ss[8]);
    br.hEnd = boost::lexical_cast<uint32_t>(ss[9]);
    br.strand = "+";
    if(br.hBeg > br.hEnd) {
        uint32_t tmp = br.hBeg;
        br.hBeg = br.hEnd;
        br.hEnd = tmp;
        br.strand = "-";
    }
    br.e = boost::lexical_cast<double>(ss[10]);
    br.score = boost::lexical_cast<float>(ss[11]);
    return br;
}

void blast_tiling(const BlastRecords& brs, string& qId, ofstream& fho) {
    fho << qId << endl;
    cout << qId << endl;
}

int main( int argc, char* argv[] ) {
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

    ifstream fhi(fi.c_str());
    if(!fhi.is_open()) { cout << format("cannot read: %s\n") % fi; return false; }
    ofstream fho(fo.c_str());
    if(!fho.is_open()) { cout << format("cannot write: %s\n") % fo; return false; }
    fho << "qId\tqBeg\tqEnd\tstrand\thId\thBeg\thEnd\tqLen\thLen\tpct\te\tscore" << endl;
    
    BlastRecords brs;
    string qId_p = "";
    string line;
    while(fhi.good()) {
        getline(fhi, line);
        if(line.length() == 0) continue;
        boost::erase_all(line, " ");
        
        vector<string> ss;
        boost::split(ss, line, boost::is_any_of("\t"));
        BlastRecord br = make_blast_record(ss);
        if(br.qId == qId_p) {
            brs.push_back(br);
        } else {
            if(qId_p != "") {
                blast_tiling(brs, qId_p, fho);
            }
            brs.clear();
            brs.push_back(br);
            qId_p = br.qId;
        }
    }
    if(brs.size() > 0) {
        blast_tiling(brs, qId_p, fho);
    }

    cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
    return EXIT_SUCCESS;
}
