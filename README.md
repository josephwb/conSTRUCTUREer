#### STRUCTUREify

##### Overview
Sample SNPs randomly across multiple phylip alignments to construct replicated STRUCTURE input files.

##### Compilation
To compile, type the following in a unix prompt:

	make

##### Usage

To run, type:

    ./STRUCTUREify -f <files> [-x seed] [-n nreps] [-o outprefix] [-p popmap] [-v]

where:

    - '-f <files>' is required and can make use of wildcards e.g. '-f *.phylip'
    - other args are optional:
        -x: random number seed (default = clock)
        -n: number of replicates to perform (default = 1)
        -o: prefix for output files (default = 'STRUCTUREify')
        -p: file indicating how individuals map to populations, with lines formatted as:
                   INDIVIDUAL_ID	POPULATION_ID
            turns on STRUCTURE format for output (default = phylip)
        -v: turn on verbose output (default = false)

For help, type:

	./STRUCTUREify -h
