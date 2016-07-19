// STRUCTUREify: Sample randomly across phylip alignments to construct STRUCTURE input files

// To compile:
// g++ -Wall STRUCTUREify.cpp -std=c++11 -O3 -o STRUCTUREify
// For run options, do:
// ./STRUCTUREify -h

#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <stdlib.h>

#include <chrono>
#include <ctime>
#include <iterator>
#include <algorithm>
#include <set>

using namespace std;

#include "utils.h"
#include "write_results.h"
#include "info.h"

// print out extra junk to screen. mostly for debugging
bool verbose = false;

// Function headers

void process_args(int argc, char *argv[], vector <string> & files, int & seed, int & nreps,
    string & outname, string & popmap, string & format);
map <string, string> read_alignment (string const& filename);
void sample_sites (vector < map <string, string> > & results, map <string, string> & new_alignment,
    vector < vector <int> > & sampled_sites);
map <string, string> read_pop_map (string & infile);

int main (int argc, char *argv[]) {
    
    int seed = -1;
    int nreps = 1;
    string outname = "STRUCTUREify";
    string mapfile = "";
    string format = "phylip";
    
    vector < map <string, string> > results;
    vector < vector <int> > sampled_sites; // for logging
    vector <string> files;
    map <string, string> pop_map;
    
    // process args (seed, nreps, files) here
    process_args(argc, argv, files, seed, nreps, outname, mapfile, format);
    
    if (format == "structure") {
        pop_map = read_pop_map(mapfile);
        cout << "Output will be in STRUCTURE format." << endl << endl;
    } else {
        cout << "Output will be in phylip format." << endl << endl;
    }
    
    cout << "Sampling " << files.size() << " files for " << nreps << " replicates." << endl;
    
    // initialize results object
    map <string, string> empty_alignment;
    for (int i = 0; i < nreps; i++) {
        results.push_back(empty_alignment);
    }
    
    // loop over files
    for (unsigned int i = 0; i < files.size(); i++) {
        map <string, string> align = read_alignment(files[i]);
        sample_sites(results, align, sampled_sites);
    }
    
    // write results
    if (format == "phylip") {
        write_results_phylip(results, outname);
    } else {
        write_results_structure(results, outname, pop_map);
    }
    
    // log which sites where sampled
    write_sampling_log(sampled_sites, files, outname);
    
    cout << endl << "Fin." << endl;
    return 0;
}



// *** LAS FUNCIONES *** //
void process_args(int argc, char *argv[], vector <string> & files, int & seed, int & nreps,
    string & outname, string & popmap, string & format) {
    program_info();
    if (argc == 1) {
        usage();
        exit(0);
    } else {
        for (int i = 1; i < argc; i++) {
            string temp = argv[i];
            
            if (temp == "-h" || temp == "-help") {
                cout <<
                "Description: sample randomly across phylip alignments to construct STRUCTURE input files" << endl <<
                endl;
                usage();
                exit(0);  
            } else if (temp == "-f") {
                i++;
                while (i < argc) {
                    string indfile = argv[i];
                    if (indfile[0] != '-') {
                        files.push_back(indfile);
                    } else {
                        i--;
                        break;
                    }
                    i++;
                }
                continue;
            } else if (temp == "-x") {
                i++;
                seed = stoi(argv[i]);
                continue;
            } else if (temp == "-n") {
                i++;
                nreps = stoi(argv[i]);
                continue;
            } else if (temp == "-o") {
                i++;
                outname = argv[i];
                continue;
            } else if (temp == "-p") {
                i++;
                popmap = argv[i];
                format = "structure";
                continue;
            } else if (temp == "-v") {
                verbose = true;
                continue;
            } else {
                cout <<
                "Unknown command-line argument '" << argv[i] << "' encountered." << endl << endl;
                usage();
                exit(0);
            }
            cout << endl;
        }
    }
    // get sorting/seeding out of the way here
    sort(files.begin(), files.end());
    seed_rand_generator(seed);
}

map <string, string> read_alignment (string const& filename) {
    string line;
    int ntax = 0;
    int nchar = 0;
    int counter = 0;
    bool first = true;
    map <string, string> align;
    
    ifstream ifs(filename.c_str());
    
    while (getline(ifs, line)) {
        if (!line.empty()) {
            vector <string> tokens = tokenize(line);
            if (first) {
                // grab dimensions: taxa char
                if (tokens.size() != 2) {
                    cout << "Error reading phylip dimensions: expected 2 values, got " <<
                        tokens.size() << ". Exiting." << endl;
                    exit (0);
                }
                ntax = stoi(tokens[0]);
                nchar = stoi(tokens[1]);
                if (verbose) {
                    cout << "Expecting " << ntax << " taxa and " << nchar << " characters for file '" <<
                        filename << "'." << endl;
                }
                
                first = false;
                continue;
            } else {
                // read in seqs
                string taxlab = prune_taxon_label(tokens[0]);
                string seq = tokens[1];
                if (seq.size() != (unsigned)nchar) {
                    cout << "Expecting " << nchar << " characters, got " << seq.size() <<
                        ". Exiting." << endl;
                }
                align[taxlab] = seq;
                counter++;
            }
        }
    }
    if (counter != ntax) {
        cout << "Error reading alignment: expected " << ntax << " taxa, got " <<
            counter << ". Exiting." << endl;
        exit(0);
    }
    return align;
}

void sample_sites (vector < map <string, string> > & results, map <string, string> & new_alignment,
    vector < vector <int> > & sampled_sites) {
    
    int nreps = results.size();
    bool first_file = false;
    if (results[0].size() == 0) {
        first_file = true;
    }
    int nchar = new_alignment.begin()->second.size();
    
    // sample sites positions
    vector <int> sites;
    for (int i = 0; i < nreps; i++) {
        sites.push_back(random_int_range(0, (nchar - 1)));
        if (verbose) {
            cout << "Sampled site " << sites[i]+1 << " (out of " << nchar << ") for rep #" << i+1 << endl;	
        }
    }
    
    // now sample actual sites, add to existing alignment
    map <string, string>::iterator iter;
    string tax;
    string seq;
    string new_char;
    for (int i = 0; i < nreps; i++) {
        for (iter = new_alignment.begin(); iter != new_alignment.end(); iter++) {
            tax = iter->first;
            seq = iter->second;
            new_char = seq[sites[i]];
            if (first_file) {
                // results are empty, nothing to match, so just add
                results[i][tax] = new_char;
            } else {
                if (results[i].count(tax) == 0) {
                    cout << "Error: taxon '" << tax << "' not found previously. Exiting." << endl;
                } else {
                    results[i][tax] += new_char;
                }
            }
        }
    }
    
    // keep for logging purposes
    sampled_sites.push_back(sites);
}





