#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cstdlib>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <Sequence/FastaExplicit.hpp>
#include <Sequence/SimpleSNP.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/PolyTableSlice.hpp>
#include <Sequence/SeqExceptions.hpp>

using namespace std;
using namespace Sequence;
using namespace Sequence::Alignment;
using namespace boost::filesystem;

path dirHome( getenv("home") );
path dirData( getenv("data") );
path dirPre( getenv("prefix") );
path dirMisc1 = dirData / "misc1";
path dirMisc2 = dirData / "misc2";

