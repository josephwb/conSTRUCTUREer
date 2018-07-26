#include <chrono>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>

using namespace std;

#include "utils.h"

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

vector <string> read_input_list_file (string const& filename) {
    vector <string> fnames;
    string line;
    ifstream ifs(filename.c_str());
    while (getline(ifs, line)) {
        if (!line.empty()) {
            fnames.push_back(line);
        }
    }
    return fnames;
}
