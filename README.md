[![DOI](https://zenodo.org/badge/256243360.svg)](https://zenodo.org/badge/latestdoi/256243360)

# N-Most-Diverse-Sequences

Summary of script purpose: The program reads in a .fasta or .afa file of sequences and, following an alignment of the sequences and creation of a distance matrix, can do any combination of the following options:

- Create a file of aligned sequences
- Produce a phylogenic tree of the sequences
- Create a file containing the names of the 'n' most diverse sequences (the user specifies the value of n if the default of 1   is not desired)

How to run: The script is run on the command line with the following command and optional flags:
  
     Rscript n_seq_script.R -i <input_file_name> 

Flags:
1. -h ...Help – display all flags
2. -v ...Version – display version of script
3. -i <file_name>...Input – required argument to specify input .fasta or .afa file
4. -t ...Tree – produces .jpeg phylogenic tree
5. -a ...Aligned - writes file of aligned sequences
6. -d ...Diversity - writes file with the names of the 'n' most diverse sequences
7. -n <integer_value>...Num - number of sequences desired for diversity methods. The default is 1.
8. -p <string_for_file_name>...Project - all files created will start with '<argument_user_entered>...'


Principle software used: R Version 3.6.3

Libraries and versions:

- getopt: version 1.20.3
- BiocManager: version 1.30.10
- seqinr: version 3.6-1
- ape: version 5.3
- DECIPHER: version 2.14.0

Principle inputs (and where files live):
- .fasta or .afa file

Created for publication (dates of running):
- (Not intended for publication)

Where software is used:
- Any server/desktop with R Version 3.6.3

Where published outputs live (if created for publication):
- (Not intended for publication)

Problems encountered while writing/running:
- Slight problems with outputting the most diverse sequences when the number of sequences desired is 1. Code outputs the correct sequence name followed by a number that should be ignored.

Project links:
- https://github.com/pfrender-laboratory/N-Most-Diverse-Sequences

Cite:

 - michalekb. 2020. “Michalekb/N-Most-Diverse-Sequences: N-Most-Diverse-Sequences,” April. https://doi.org/10.5281/ZENODO.3754619.
