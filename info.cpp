// program info

#include <iostream>

using namespace std;

#include "info.h"

// version information
string name = "STRUCTUREify";
double version = 0.1;
string month = "July";
int year = 2016;

void program_info() {
    cout << endl << 
    "************************************************" << endl <<
    "           STRUCTUREify version " << version      << endl <<
    "                Joseph W. Brown"                  << endl <<
    "             University of Michigan"              << endl <<
    "         Complaints: josephwb@umich.edu"          << endl <<
    "                  " << month <<", " << year       << endl << 
    "************************************************" << endl << endl;
}

void usage () {
    cout << "Usage:" << endl << endl <<
        "./STRUCTUREify -f <files> [-x seed] [-n nreps] [-o outprefix] [-p popmap] [-v]" << endl << endl <<
        "where:" << endl << 
        "- '-f <files>' is required and can make use of wildcards e.g. '-f *.phylip'" << endl <<
        "- '-l <list>' is a file containing names of input files, 1 per line" << endl <<
        "    - use this option if you get a 'Argument list too long' error message" << endl <<
        "- other args are optional:" << endl <<
        "    -x: random number seed (default = clock)" << endl <<
        "    -n: number of replicates to perform (default = 1)" << endl <<
        "    -o: prefix for output files (default = 'STRUCTUREify')" << endl <<
        "    -p: file indicating how individuals map to populations, with lines formatted as:" << endl <<
        "               INDIVIDUAL_ID	POPULATION_ID" << endl <<
        "        turns on STRUCTURE format for output (default = phylip)" << endl <<
        "    -v: turn on verbose output (default = false)" << endl << endl;
}
