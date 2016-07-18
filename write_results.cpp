#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <ctime>
#include <set>

using namespace std;

#include "write_results.h"
#include "utils.h"

extern bool verbose;

// this is how nucleotides are mapped to integer alleles.
// NOTE: indel is currently a fifth state (i.e. not missing)
// TODO: make this malleable
map <char, int> nucToInt = {
    {'A', 1},
    {'T', 2},
    {'G', 3},
    {'C', 4},
    {'N', -9},
    {'-', 5}
};

void write_results_phylip (vector < map <string, string> > & results, string & outname) {
    int nreps = results.size();
    
    map <string, string>::iterator iter;
    for (int i = 0; i < nreps; i++) {
        
        ofstream outFile;
        string ofname = outname;
        if (nreps > 1) {
            ofname += "_rep" + to_string(i+1) + ".phy"; 
        } else {
            ofname += ".phy";
        }
        
        outFile.open(ofname.c_str());
        
        if (verbose) {
            cout << endl << "Results for sampling replicate #" << i+1 << ":" << endl;
            cout << results[i].size() << " " << results[i].begin()->second.size() << endl;
        }
        outFile << results[i].size() << " " << results[i].begin()->second.size() << endl;
        for (iter = results[i].begin(); iter != results[i].end(); iter++) {
            string lab = iter->first;
            int diff = 10 - lab.size();
            lab += string(diff, ' ');
            if (verbose) {
                cout << lab << "  " << iter->second << endl;
            }
            outFile << lab << "  " << iter->second << endl;
        }
        outFile.close();
    }
    if (nreps == 1) {
        cout << endl << "Wrote results to file '" << outname << ".phy'." << endl;
    } else {
        cout << endl << "Wrote " << nreps << " result files '" << outname << "_repN.phy" << "'." << endl;
    }
}


void write_results_structure(vector < map <string, string> > & results, string & outname,
    map <string, string> & pop_map) {
    int nreps = results.size();
    
    map <string, string>::iterator iter;
    for (int i = 0; i < nreps; i++) {
        
        ofstream outFile;
        string ofname = outname;
        if (nreps > 1) {
            ofname += "_rep" + to_string(i+1) + ".str"; 
        } else {
            ofname += ".str";
        }
        
        outFile.open(ofname.c_str());
        
        if (verbose) {
            cout << endl << "Results for sampling replicate #" << i+1 << ":" << endl;
        }
        for (iter = results[i].begin(); iter != results[i].end(); iter++) {
            string lab = iter->first;
            string seq = iter->second;
            int diff = 10 - lab.size();
            lab += string(diff, ' ');
            if (verbose) {
                cout << lab << "	" << pop_map[iter->first] << "	";
            }
            outFile << lab << "	" << pop_map[iter->first] << "	";
            for (unsigned int j = 0; j < seq.size(); j++) {
                outFile << to_string(nucToInt[seq[j]]);
                if (j < (seq.size() - 1)) {
                     outFile << " ";
                }
                if (verbose) {
                    cout << to_string(nucToInt[seq[j]]);
					if (j < (seq.size() - 1)) {
						 cout << " ";
					}
				}
            }
            if (verbose) {
                cout << endl;
            }
            outFile << endl;
        }
        outFile.close();
    }
    if (nreps == 1) {
        cout << endl << "Wrote results to file '" << outname << ".str'." << endl;
    } else {
        cout << endl << "Wrote " << nreps << " result files '" << outname << "_repN.str" << "'." << endl;
    }
}

void write_sampling_log (vector < vector <int> > const& sampled_sites, vector <string> const& files,
    string & outname) {
    int nreps = sampled_sites[0].size();
    int nfiles = files.size();
    
    string ofname = outname + "_sampling.log";
    
    ofstream outFile;
    outFile.open(ofname.c_str(), ios::out | ios::app);
    
    time_t t = time(NULL);
    char* charTime = ctime(&t);
    outFile << "# " << charTime << endl;
    
    cout << endl << "Sampled sites (file	rep1_site rep2site etc.):" << endl << endl;
    
    for (int i = 0; i < nfiles; i++) {
        cout << files[i] << "	";
        outFile << files[i] << "	";
        for (int j = 0; j < nreps; j++) {
            cout << sampled_sites[i][j] + 1;
            outFile << sampled_sites[i][j] + 1;
            if (j < (nreps - 1)) {
                cout << " ";
                outFile << " ";
            }
        }
        cout << endl;
        outFile << endl;
    }
    outFile << endl;
    outFile.close();
    
    cout << endl << "Sampling info written to log file '" << ofname << "'." << endl;
}

map <string, string> read_pop_map (string & infile) {
    map <string, string> popMap;
    string line;
    set <string> uniquePops;
    ifstream ifs(infile.c_str());
    cout << endl << "Reading in population mapping file '" << infile << "'." << endl;
    while (getline(ifs, line)) {
        if (!line.empty()) {
            vector <string> tokens = tokenize(line);
            popMap[tokens[0]] = tokens[1];
            uniquePops.insert(tokens[1]);
        }
    }
    cout << "Read in mapping of " << popMap.size() << " individuals to " << uniquePops.size() <<
        " populations:" << endl << endl;
    
    cout << "INDIVIDUAL  POP" << endl;
    for (map <string, string>::iterator iter = popMap.begin(); iter != popMap.end(); iter++) {
        string lab = iter->first;
		int diff = 10 - lab.size();
		lab += string(diff, ' ');
        cout << lab << "  " << iter->second << endl;
    }
    cout << endl;
    
    return popMap;
}
