#ifndef _WRITE_RES_
#define _WRITE_RES_

void write_results_phylip (vector < map <string, string> > & results, string & outname);
void write_results_structure(vector < map <string, string> > & results, string & outname,
    map <string, string> & pop_map);
void write_sampling_log (vector < vector <int> > const& sampled_sites, vector <string> const& files,
    string & outname);

#endif /* _WRITE_RES_ */
