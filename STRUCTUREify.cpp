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

// this is how nucleotides are mapped to integer alleles.
// NOTE: indel is currently a fifth state (i.e. not missing)
map <char, int> nucToInt = {
    {'A', 1},
    {'T', 2},
    {'G', 3},
    {'C', 4},
    {'N', -9},
    {'-', 5}
};

// version information
string name = "STRUCTUREify";
double version = 0.1;
string month = "July";
int year = 2016;

// print out extra junk to screen. mostly for debugging
bool verbose = false;

// Function headers
void program_info();
void usage();
void process_args(int argc, char *argv[], vector <string> & files, int & seed, int & nreps,
    string & outname, string & popmap, string & format);
void seed_rand_generator (int const& seed);
int random_int_range (int const& min, int const& max);
vector <string> tokenize (string const& input);
string prune_taxon_label (string & label);
map <string, string> read_alignment (string const& filename);
void sample_sites (vector < map <string, string> > & results, map <string, string> & new_alignment,
    vector < vector <int> > & sampled_sites);
void write_results_phylip (vector < map <string, string> > & results, string & outname);
void write_results_structure(vector < map <string, string> > & results, string & outname,
    map <string, string> & pop_map);
void write_sampling_log (vector < vector <int> > const& sampled_sites, vector <string> const& files,
    string & outname);
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
        "- other args are optional:" << endl <<
        "    -x: random number seed (default = clock)" << endl <<
        "    -n: number of replicates to perform (default = 1)" << endl <<
        "    -o: prefix for output files (default = 'STRUCTUREify')" << endl <<
        "    -p: file indicating how individuals map to populations, with lines formatted as:" << endl <<
        "               INDIVIDUAL_ID	POPULATION_ID" << endl <<
        "        turns on STRUCTURE format for output (default = phylip)" << endl <<
        "    -v: turn on verbose output (default = false)" << endl << endl;
        
}

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

void seed_rand_generator (int const& seed) {
    if (seed == -1) {
        srand(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    } else {
        srand(seed);
    }
}

int random_int_range (int const& min, int const& max) {
    return min + (rand() % (int)(max - min + 1));
}

vector <string> tokenize (string const& input) {
    vector <string> tokens;
    string temp;
    istringstream str(input);
    while (str >> temp) {
        tokens.push_back(temp);
    }
    return tokens;
}

string prune_taxon_label (string & label) {
    // prune everything after (and including) first '_'
    string::size_type pos = label.find_first_of('_');
    //cout << "Original label: " << label;
    if (string::npos != pos) {
        label.erase(label.begin() + pos, label.end());
    }
    //cout << "; New label: " << label << endl;
    return label;
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



