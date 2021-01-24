# search_TPS - A script to search for monoTPS, diTPS and sesquiTPS sequences in a set of input proteins.

search_TPS is a computational tool for the detection and classification of TPS sequences of the monoTPS, diTPS and sesquiTPS classes using class-specific profile PFAMs and HMMs.

The program receives as input a file containing the sequences of interest in FASTA format, a directory composed of the PFAM models of terpene synthase, a directory composed of class-specific profile HMMs, and a tabular file containing the cut-off score values for each class-specific profile HMM. To identify low confidence TPS sequences, the program uses a standard score value (score = 100), but this value can be set by the user.

From the cutoff score values of the class-specific profile HMMs, the program identifies highly confidence TPS sequences contained in the input file. To select low confidence TPS sequences, the program uses a standard score value (score = 100), but this value can be set by the user.

##   Instalation

search_TPS does not need to be installed. The user should only download the search_TPS.pl file

## Requirements

search_TPS only requires the program hmmsearch from HMMER3 package (http://hmmer.org/) to perform similarity searches of the profile HMMs against the input sequences. The hmmsearch  program must be located in a directory listed in the PATH of the operating system.

## Usage
```
search_TPS.pl -d <class-scpecific_phmms_dir> -i <input_fasta_file> -t <table_file> -s <search_option> -p <pfams_dir> -o <output_dir>
```
### Example:
```
search_TPS.pl -d class-specific_hmms -i sequences.fasta -t score_tables_dir/all.scores -s 1 -p PFAMs_dir -o results_dir
```

### Mandatory parameters:
```
-d <class-scpecific_phmms_dir>  : Directory containing specific profile HMMs for monoTP, diTP and sesquiTP.
-i <input_fasta_file>           : Fasta file containing the sequences to be searched.
-p <pfams_dir>                  : Directory containing the PFAMs models.
-s <search_options>             : Group of TPS to be searched:
                                        1 - all (monoTP, diTP and sesquiTP)
                                        2 - monoTP
                                        3 - diTP
                                        4 - sesquiTP
-t <table_file>                 : File containing the cutoff scores to be used to select the results of each HMM profile.
```

### Optional parameters:
```
-cpu <num>                      : Number of threads to be used by hmmsearch.
-h|help                         : Show this help message.
-l <low_confidence_score>       : Low confidence score to be considered (default = 100).
-o <dir>                        : Output directory (default = output_dir).
-v|version                      : Version.
```

## Contact

To report bugs, to ask for help and to give any feedback, please contact Liliane S. Oliveira (liliane.sntn@gmail.com).
