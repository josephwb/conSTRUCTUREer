#ifndef _UTILS_
#define _UTILS_

void seed_rand_generator (int const& seed);
int random_int_range (int const& min, int const& max);
vector <string> tokenize (string const& input);
string prune_taxon_label (string & label);

#endif /* _UTILS_ */
