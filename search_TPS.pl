#!/usr/bin/perl

# search_TPS.pl v. 1.0 (2021-01-23)

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

#variables
my $output = 'output_dir';
my $version = "1.0";
my $last_update = "2021-01-23";
my $help;
my $hmms_dir;
my $table;
my $pfams_dir;
my $low_score = 100;
my $hmmsearch_tabular_file;
my $transeq_file;
my $option;
my $cpu;
my $input_file;
my $file_log = "file.log";
my $version_op;
my $help_print = "##### Search TPS - version $version ($last_update) - L. Oliveira#####


Usage:
search_TPS.pl -d <class-scpecific_phmms_dir> -i <input_fasta_file> -t <table_file> -s <search_option> -p <pfams_dir> -o <output_dir> 

Mandatory parameters:
-d <class-scpecific_phmms_dir>	: Directory containing specific profile HMMs for monoTP, diTP and sesquiTP.
-i <input_fasta_file> 		: Fasta file containing the sequences to be searched.
-p <pfams_dir>			: Directory containing the PFAMs models. 
-s <search_options>		: Group of TPS to be searched:
					1 - all (monoTP, diTP and sesquiTP)
					2 - monoTP
					3 - diTP
					4 - sesquiTP
-t <table_file>         	: File containing the cutoff scores to be used to select the results of each HMM profile.

OPTIONAL PARAMETERS:
-cpu <num>		  	: Number of threads to be used by hmmsearch.
-h|help           		: Show this help message.
-l <low_confidence_score>       : Low confidence score to be considered (default = 100).
-o <dir>            	  	: Output directory (default = output_dir).
-v|version        		: Version.
 \n";

my $optret = GetOptions ("d=s"			=> \$hmms_dir,
                         "i=s"  	        => \$input_file,
                         "s=i"                  => \$option,
                         "o=s"		        => \$output,
			 "p=s"			=> \$pfams_dir,
			 "l=f"			=> \$low_score,
                         "h|help"               => \$help,
			 "t=s"			=> \$table,
			 "cpu=i"		=> \$cpu,
			 "v|version"            => \$version_op);

if($help){
    print $help_print;
    die "\n";    
}

if($version_op){
    die "Version $version.\nLast update: $last_update\n";
}

if(!$input_file){
    die "ERROR: Missing mandatory argument -i.\n$help_print\n";
}

if(!$hmms_dir){
    die "ERROR: Missing mandatory argument -d.\n$help_print\n";
}

if(!$option){
    die "ERROR: Missing mandatory argument -s
    .\n$help_print\n";
}

if(!$table){
    die "ERROR: Missing mandatory argument -t.\n$help_print\n";
}

if(!$cpu){
    $cpu = `grep CPU -c /proc/cpuinfo`;
    $cpu = $cpu/2;
}

# Creates the output directory
$output = output_dir_name($output);
my $original_dir = $output;
system "mkdir $output";
$file_log = $output."/".$file_log;

# Creates the log file
my $log = "$output/error.log";
open(my $log_file_handle, ">$log") or die "ERROR: Could not create log file $log : $!\n";

my $transeq_name;
open(my $fl, ">$file_log") or die "ERROR: Could not create log file $file_log : $!\n";
print $fl "Input fasta file: $input_file\n";
print $fl "Profile HMMs directory: $hmms_dir\n";

# Stores PFAM models in memory
opendir(DIR, "$pfams_dir");
my @pfams = readdir(DIR);
closedir(DIR);

my @models = ();
my @monoTP = ();
my @diTP = ();
my @sesquiTP = ();
my $search_type;

# Stores the class-specific profile HMMs in memory
if($option == 1){
    $search_type = "MonoTP, diTP and sesquiTP";
    my $aux = `ls $hmms_dir/*monoTP.hmm`;
    @monoTP = split("\n", $aux);
    $aux = `ls $hmms_dir/*diTP.hmm`;
    @diTP = split("\n", $aux);
    $aux = `ls $hmms_dir/*sesquiTP.hmm`;
    @sesquiTP = split("\n", $aux);
}
elsif($option == 2){
    $search_type = "MonoTP";
    my $aux = `ls $hmms_dir/*monoTP.hmm`;
    @models = split("\n", $aux);
}
elsif($option == 3){
    $search_type = "diTP";
    my $aux = `ls $hmms_dir/*diTP.hmm`;
    @models = split("\n", $aux);
}
elsif($option == 4){
    $search_type = "sesquiTP";
    my $aux = `ls $hmms_dir/*sesquiTP.hmm`;
    @models = split("\n", $aux);
}
else{
    die "Invalid option! Select 0, 1, 2 or 3!\n";
}
print $fl "Searching for: $search_type.\n";
print $fl "Scores table: $table\n";

my %scores = ();
open(FILE, $table);
while(<FILE>){
    chomp($_);
    if($_ =~ /Profile HMM/){
	next;
    }
    my @aux = split("\t", $_);
    $scores{$aux[0]} = $aux[1];
}
close(FILE);

# Verifies the composition of the fasta file (nucleotide or amino acid)
my $type = verifyFastaFileComposition($input_file);
if($type == 1){ # Nucleotide
    print $fl "Type of the input fasta file: Nucleotide Fasta\n";
    # Runs transeq program to translate the nucleotide sequences into proteins
    $transeq_name = runTranseq($input_file, $output);
}
else{# Amino acid
    print $fl "Type of the input fasta file: Protein Fasta\n";
    $transeq_name = $input_file;
}

my %seqs = ();
# Translates the input nucleotide sequences to protein using the 6 frames
open(FILE, $transeq_name);
my $seq = "";
my $name = "";
while(<FILE>){
    chomp($_);
    if($_ =~ /^>/){
	if($seq ne ""){
	    $seqs{$name} = $seq;
	}
	$seq = "";
	my @aux = split(" ",$_);
	$name = $aux[0];
	$name =~ s/>//g;
    }
    else{
	$seq .= $_;
    }
}
$seqs{$name} = $seq;
close(FILE);

# Runs similarity searches of the PFAM models against the input file and selects the positive sequences
print STDERR "Selecting sequences using PFAM models...\n";
runPFAMsSimilaritySearch($pfams_dir, "PFAM", $transeq_name);

my $selected = $output."/PFAM_dir/selected.fasta";
my $pfam_selected = $output."/PFAM_dir/final_results.csv";
open(OUT, ">$selected");
open(FILE, $pfam_selected);
while(<FILE>){
    chomp($_);
    print OUT ">$_\n$seqs{$_}\n";
}
close(FILE);
close(OUT);


if($option == 1){ # Runs semilarities searches of the class-specific profile HMMs for all classes against the select sequences
    print STDERR "Searching for monoTP sequences...\n";
    runSimilaritySearch(\@monoTP, "monoTP", $selected);
    print STDERR "Searching for diTP sequences...\n";
    runSimilaritySearch(\@diTP, "diTP", $selected);
    print STDERR "Searching for sesquiTP sequences...\n";
    runSimilaritySearch(\@sesquiTP, "sesquiTP", $selected);
    # generate final results
    removeRedundantSequences($output);
}
else{ # Runs semilarities searches of the class-specific profile HMMs for a defined class against the select sequences
    my $out_dir = $output."/".$search_type."_dir";
    system "mkdir $out_dir";
    for my $hmm (@models){
    	if($hmm eq "." or $hmm eq ".."){
	   next;
    	}
    	my $pos = rindex($hmm, "/");
    	my $name = substr $hmm, ($pos+1);
    	print STDERR "Profile HMM $name...\n";
    	my $full = "$out_dir/".$name."_hmmsearch.txt";    
    	open(FILE, ">$full");
    	close(FILE);
    	my $low_tab = "$out_dir/".$name."_hmmsearch_low.tab";
    	open(FILE, ">$low_tab");
    	close(FILE);
	my $high_tab = "$out_dir/".$name."_hmmsearch_high.tab";
        open(FILE, ">$high_tab");
        close(FILE);
    	runHmmsearch($hmm, $selected, $full, $low_tab, $high_tab, $output);
    }
    print STDERR "Analysing the results...\n";
    generateFinalResultsOneGroup($out_dir);
}

print "Done.\n";
print $fl "Done.\n";
close($fl);

exit;

################################################################
### Sobroutines
################################################################

# Runs similarity searches using PFAM models
sub runPFAMsSimilaritySearch{
    my $dir = shift;
    my $name = shift;
    my $file = shift;
    opendir(DIR, "$dir");
    my @pfams = readdir(DIR);
    closedir(DIR);
    my $out_dir = $output."/".$name."_dir";
    system "mkdir $out_dir";
    for my $hmm (@pfams){
        if($hmm eq "." or $hmm eq ".."){
           next;
        }
        my $pos = rindex($hmm, "/");
        my $name = substr $hmm, ($pos+1);
        print STDERR "Profile HMM $name...\n";
        my $full = "$out_dir/".$name."_hmmsearch.txt";
        open(FILE, ">$full");
        close(FILE);
        my $tab = "$out_dir/".$name."_hmmsearch.tab";
        open(FILE, ">$tab");
        close(FILE);
	my $model = $dir."/$hmm";
	my $command = "hmmsearch --tblout $tab -o $full --cpu $cpu $model $file 2>> $log";
        system $command;	
    }
    print STDERR "Analysing the results...\n";
    generateFinalResultsPFAMs($out_dir);
}

# Runs similarity searches using class-specific profile HMMs
sub runSimilaritySearch{
    my $array = shift;
    my $name = shift;
    my $file = shift;
    my $dir = $output."/".$name."_dir";
    system "mkdir $dir";
    my @models = @{$array};
    for my $hmm (@models){
        if($hmm eq "." or $hmm eq ".."){
           next;
        }
        my $pos = rindex($hmm, "/");
        my $name = substr $hmm, ($pos+1);
        print STDERR "Profile HMM $name...\n";
        my $full = "$dir/".$name."_hmmsearch.txt";
        open(FILE, ">$full");
        close(FILE);
        my $low_tab = "$dir/".$name."_hmmsearch_low.tab";
        open(FILE, ">$low_tab");
        close(FILE);
	my $high_tab = "$dir/".$name."_hmmsearch_high.tab";
        open(FILE, ">$high_tab");
        close(FILE);
        runHmmsearch($hmm, $file, $full, $low_tab, $high_tab, $output);
    }
    print STDERR "Analysing the results...\n";
    generateFinalResults($dir);
}

# Runs similarity searches
sub runHmmsearch{
    my $hmm = shift;
    my $fasta = shift;
    my $full = shift;
    my $low_tab = shift;
    my $high_tab = shift;
    my $out = shift;
    my $dir = $out."/hmms";
    system "mkdir $dir";
    my $results = $out."/results";
    system "mkdir $results";
    $/ = "//\n";
    open(DATA, "<$hmm");
    my $count = 0;
    foreach my $unique_hmm (<DATA>){
    	my @hmm_lines = split("\n", $unique_hmm);
	my $aux_n;
	my $achei = 0;
        foreach my $line (@hmm_lines){
	    if($line =~ m/(NAME\d*)/g){
	        $aux_n = $line;
		$aux_n =~ s/NAME//;
		$aux_n =~ s/\s+//g;
		last;
	    }
	}
        my $file = $dir."/".$aux_n.".hmm";
	open(HMMFILE, ">$file") or die "ERROR: Couldn't create file $output/$file : $!\n";
        print HMMFILE $unique_hmm;
        close (HMMFILE);
    }
    close(DATA);
    $/ = "\n";
    opendir(DIR, "$dir");
    my @new_models = readdir(DIR);
    closedir(DIR);
    open(FULL, ">$full");
    open(LOW, ">$low_tab");
    open(HIGH, ">$high_tab");
    for my $model (@new_models){
	if($model eq "." or $model eq ".."){
	    next;
	}
	my $aux_n = $model;
	$aux_n =~ s/.hmm//;
        my $file = $dir."/".$model;	    
	my $file_full = $results."/".$aux_n."_full";
        my $file_tab = $results."/".$aux_n."_tab";
	my $score = 1;
	if(defined $scores{$aux_n}){
	    $score = $scores{$aux_n};
	}
	my $command = "hmmsearch --tblout $file_tab -o $file_full --cpu $cpu $file $fasta 2>> $log";
	system $command; 
	open(FILE, $file_full);
        while(<FILE>){
            chomp($_);
            if($_ =~ /#/){
            	if($count == 0){
                    print FULL "$_\n";
                }
            }
            else{
                print FULL "$_\n";
            }
        }
	close(FILE);
        open(FILE, $file_tab);
        while(<FILE>){
            chomp($_);
            if($_ =~ /#/){
            	if($count == 0 and ($_ =~ /full sequence/ or $_ =~ /target name/)){
                    print HIGH "$_\n";
		    print LOW "$_\n";
                }
            }
            else{
		my @aux = split(" ", $_);
		if($aux[5] >= $score){
		    print HIGH "$_\n";
		    print LOW "$_\n";
		}
		elsif($aux[5] >= $low_score and $aux[5] < $score){
                    print LOW "$_\n";
		}
            }
        }
        close(FILE);
        ++$count;
    }
    close(FULL);
    close(LOW);
    close(HIGH);
    system "rm -rf $results $dir";
}

# Checks if a directory with the given output name already exists. If so, a numeric suffix is added to the name.
sub output_dir_name {
    my $output_dir_name = shift;
    my $count = 2;
    my $flag = 0;
    print STDERR "Creating output directory $output_dir_name\n";
    if (-d $output_dir_name) {
        $flag = 1;
        while (-d "$output_dir_name\_$count") {
            $count++;
        }
        $output_dir_name = "$output_dir_name\_$count";
    }
    print STDERR "\nOutput directory already exists, saving results to $output_dir_name instead.\n\n" if $flag;
    return ($output_dir_name);
}

# Identifies whether the fasta file contains nucleotide or protein sequences.
sub verifyFastaFileComposition{
    my $file = shift;
    my $count = 0;
    my $num_lines = 0;
    my $type = 1; # 1 - nucleotide; 2 - protein
    open(my $vf, "$file");
    while(<$vf>){
    	chomp($_);
	if($_ =~ /^>/){}
        else{
	    if ((index($_, 'L') != -1) or (index($_, 'l') != -1)) {
	    	++$count;
	    }
	    if ((index($_, 'I') != -1) or (index($_, 'i') != -1)) {
            	++$count;
            }
	    if ((index($_, 'P') != -1) or (index($_, 'p') != -1)) {
            	++$count;
            }
	    if ((index($_, 'F') != -1) or (index($_, 'f') != -1)) {
            	++$count;
            }
	    if ((index($_, 'Q') != -1) or (index($_, 'q') != -1)) {
            	++$count;
            }
	    if ((index($_, 'E') != -1) or (index($_, 'e') != -1)) {
            	++$count;
            }
	    if($count > 0){
	    	$type = 2;
	    	last;
	    }
	    else{
		++$num_lines;
		if($num_lines > 20){
		    last;
		}
	    }
        }
    }
       
    close($file);
    return $type;
}

# Runs transeq program.
sub runTranseq{
    my $file = shift;
    my $dir = shift;
    my $string;
    my $aux;    
    if(defined $dir){
    	$string = $dir."/".$file."*_transeq*";
	$aux = `ls $string 2> /dev/null`;
	if($aux eq ""){
	    $string = $file."*_transeq*";
	    $aux = `ls $string 2> /dev/null`;
	}
    }
    else{
	$string = $file."*_transeq*";
	$aux = `ls $string 2> /dev/null`;	
    }

    if($aux eq ""){	
	$aux = undef;
    }
    else{
	my @aux_t = split("\n", $aux);
	$aux = $aux_t[0];
    }
    my $transeq;
    my $auxiliar;
    print $fl "Step: translation of nucleic acid sequences\n";
    if(defined $aux){ #Verify if transeq file exists
	print $fl "Transeq protein file found. Skipping translation of nucleic acid sequences\n";
	print "Transeq protein file found. Skipping translation of nucleic acid sequences\n";
	$transeq = $aux;
    }
    else{ #run the transeq program:
     	$transeq = $file."_transeq.fasta";
	print $fl "Performing conceptual translation of $file and creating $transeq (Parameters: frame = 6) \n";
        print "Performing conceptual translation of $file and creating $transeq... \n";
        system "transeq -frame 6 $file $transeq";
        print "Done.\n";
	print $fl "Done\n";
    }
    return $transeq;
}

# Generates final results for PFAM models: stores the sequences presenting similarity with at least two models in a new file.
sub generateFinalResultsPFAMs{
    my $dir = shift;
    my %all = ();
    my $out = $dir."/final_results.csv";
    my $aux = `ls $dir/*tab`;
    my @files = split("\n", $aux);
    for my $file (@files){
        if($file eq "." or $file eq ".."){
            next;
        }
        open(FILE, $file);
        while(<FILE>){
            chomp($_);
            if($_ =~ /^#/){
                next;
            }
            my @line = split(" ", $_);
            if(defined $all{$line[0]}){
                $all{$line[0]} += 1;
            }
            else{
                $all{$line[0]} = 1;
            }
        }
        close(FILE);
    }
    open(OUT, ">$out");
    foreach my $key (sort keys %all){
        if($all{$key} > 1){
            print OUT "$key\n";
        }
    }
    close(OUT);
}	    

# Generates final results for class-specific profile HMMs: stores the sequences 
# presenting similarity with at least two models in a new file
sub generateFinalResults{
    my $dir = shift;
    my %low = ();
    my %high = ();
    my %high_scores = ();
    my %low_scores = ();
    my $low_file = $dir."/low_confidence_all_results.csv";
    my $high_file = $dir."/high_confidence_all_results.csv";
    my $aux = `ls $dir/*_high.tab`;
    my @files = split("\n", $aux);
    for my $file (@files){
        if($file eq "." or $file eq ".."){
            next;
        }
        open(FILE, $file);
        while(<FILE>){
            chomp($_);
            if($_ =~ /^#/){
                next;
            }
            my @line = split(" ", $_);
	    if(defined $high_scores{$line[0]}){
		if($high_scores{$line[0]} < $line[5]){
	    	    $high_scores{$line[0]} = $line[5];
		}
	    }
	    else{
	    	$high_scores{$line[0]} = $line[5];
	    }
            if(defined $high{$line[0]}){
                $high{$line[0]} += 1;
            }
            else{
                $high{$line[0]} = 1;
            }
        }
        close(FILE);
    }
    $aux = `ls $dir/*_low.tab`;
    @files = split("\n", $aux);
    for my $file (@files){
        if($file eq "." or $file eq ".."){
            next;
        }
        open(FILE, $file);
        while(<FILE>){
            chomp($_);
            if($_ =~ /^#/){
                next;
            }
            my @line = split(" ", $_);
	    if(defined $low_scores{$line[0]}){
		if($low_scores{$line[0]} < $line[5]){
                    $low_scores{$line[0]} = $line[5];
		}
            }
            else{
                $low_scores{$line[0]} = $line[5];
            }
            if(defined $low{$line[0]}){
                $low{$line[0]} += 1;
            }
            else{
                $low{$line[0]} = 1;
            }
        }
        close(FILE);
    }
    open(OUT, ">$high_file");
    my %selected_high = ();
    print OUT "Sequence\tHighest score\n";
    foreach my $key (sort keys %high){
        if($high{$key} > 1){
            $selected_high{$key} = 1;
            print OUT "$key\t$high_scores{$key}\n";
        }
	elsif($high{$key} == 1 and $low{$key} >= 1){
	    $selected_high{$key} = 1;
            print OUT "$key\t$high_scores{$key}\n";
	}
    }
    close(OUT);

    open(OUT, ">$low_file");
    print OUT "Sequence\tHighest score\n";
    foreach my $key (sort keys %low){
        if($low{$key} > 1 and !defined $selected_high{$key}){
            print OUT "$key\t$low_scores{$key}\n";
        }
    }
    close(OUT);

}

# Classifies sequences presenting similarity with more than one class according their scores.
sub removeRedundantSequences{
    my $dir =  shift;
    my %mono_high = ();
    my $file = $dir."/monoTP_dir/high_confidence_all_results.csv";
    open(FILE, $file);
    while(<FILE>){
        chomp($_);
	if($_ =~ /Highest score/){
            next;
        }
        my @aux = split("\t", $_);
        $mono_high{$aux[0]} = $aux[1];
    }
    close(FILE);

    $file = $dir."/monoTP_dir/low_confidence_all_results.csv";
    my %mono_low = ();
    open(FILE, $file);
    while(<FILE>){
        chomp($_);
	if($_ =~ /Highest score/){
	    next;
	}
        my @aux = split("\t", $_);
        $mono_low{$aux[0]} = $aux[1];
    }
    close(FILE);

    $file = $dir."/diTP_dir/high_confidence_all_results.csv";    
    my %di_high = ();
    open(FILE, $file);
    while(<FILE>){
        chomp($_);
	if($_ =~ /Highest score/){
            next;
        }
        my @aux = split("\t", $_);
        $di_high{$aux[0]} = $aux[1];
    }
    close(FILE);

    $file = $dir."/diTP_dir/low_confidence_all_results.csv";
    my %di_low = ();
    open(FILE, $file);
    while(<FILE>){
        chomp($_);
        if($_ =~ /Highest score/){
            next;
        }
        my @aux = split("\t", $_);
        $di_low{$aux[0]} = $aux[1];
    }
    close(FILE);
    
    $file = $dir."/sesquiTP_dir/high_confidence_all_results.csv";
    my %sesqui_high = ();
    open(FILE, $file);
    while(<FILE>){
        chomp($_);
        if($_ =~ /Highest score/){
            next;
        }
        my @aux = split("\t", $_);
        $sesqui_high{$aux[0]} = $aux[1];
    }
    close(FILE);

    $file = $dir."/sesquiTP_dir/low_confidence_all_results.csv";
    my %sesqui_low = ();
    open(FILE, $file);
    while(<FILE>){
        chomp($_);
        if($_ =~ /Highest score/){
            next;
        }
        my @aux = split("\t", $_);
        $sesqui_low{$aux[0]} = $aux[1];
    }
    close(FILE);
    
    my $out_file = $dir."/monoTP_dir/high_confidence_final_results.csv";
    open(OUT, ">$out_file");
    print OUT "Sequence\tHighest score\n";
    my %present = ();
    foreach my $key (sort keys %mono_high){
	my $found = 0;
	if(defined $di_high{$key} and $di_high{$key} > $mono_high{$key}){
	    next;
	}
	if(defined $sesqui_high{$key} and $sesqui_high{$key} > $mono_high{$key}){
            next;
        }
	$present{$key} = $mono_high{$key};
	print OUT $key."\t".$mono_high{$key}."\n";
    }
    close(OUT);

    $out_file = $dir."/monoTP_dir/low_confidence_final_results.csv";
    open(OUT, ">$out_file");
    print OUT "Sequence\tHighest score\n";
    foreach my $key (sort keys %mono_low){
        my $found = 0;
	if(defined $present{$key}){
	    next;
	}
        if(defined $di_high{$key} and $di_high{$key} > $mono_low{$key}){
            next;
        }
        if(defined $di_low{$key} and $di_low{$key} > $mono_low{$key}){
            next;
        }
        if(defined $sesqui_high{$key} and $sesqui_high{$key} > $mono_low{$key}){
            next;
        }
        if(defined $sesqui_low{$key} and $sesqui_low{$key} > $mono_low{$key}){
            next;
        }
        $present{$key} = $mono_low{$key};
	print OUT $key."\t".$mono_low{$key}."\n";
    }
    close(OUT);

    $out_file = $dir."/diTP_dir/high_confidence_final_results.csv";
    open(OUT, ">$out_file");
    print OUT "Sequence\tHighest score\n";
    foreach my $key (sort keys %di_high){
        my $found = 0;
        if(defined $present{$key}){
            next;
        }
        if(defined $sesqui_high{$key} and $sesqui_high{$key} > $di_high{$key}){
            next;
        }
        $present{$key} = $di_high{$key};
	print OUT $key."\t".$di_high{$key}."\n";
    }
    close(OUT);

    $out_file = $dir."/diTP_dir/low_confidence_final_results.csv";
    open(OUT, ">$out_file");
    print OUT "Sequence\tHighest score\n";
    foreach my $key (sort keys %di_low){
        my $found = 0;
        if(defined $present{$key}){
            next;
        }
        if(defined $sesqui_high{$key} and $sesqui_high{$key} > $di_low{$key}){
            next;
        }
        if(defined $sesqui_low{$key} and $sesqui_low{$key} > $di_low{$key}){
            next;
        }
        $present{$key} = $di_low{$key};
	print OUT $key."\t".$di_low{$key}."\n";
    }
    close(OUT);

    $out_file = $dir."/sesquiTP_dir/high_confidence_final_results.csv";
    open(OUT, ">$out_file");
    print OUT "Sequence\tHighest score\n";
    foreach my $key (sort keys %sesqui_high){
        my $found = 0;
        if(defined $present{$key}){
            next;
        }
        $present{$key} = $sesqui_high{$key};
	print OUT $key."\t".$sesqui_high{$key}."\n";
    }
    close(OUT);

    $out_file = $dir."/sesquiTP_dir/low_confidence_final_results.csv";
    open(OUT, ">$out_file");
    print OUT "Sequence\tHighest score\n";
    foreach my $key (sort keys %sesqui_low){
        my $found = 0;
        if(defined $present{$key}){
            next;
        }
        $present{$key} = $sesqui_low{$key};
	print OUT $key."\t".$sesqui_low{$key}."\n";
    }
    close(OUT);
}

# Generates final results for class-specific profile HMMs (when only one class is searched): stores the sequences
# presenting similarity with at least two models in a new file
sub generateFinalResultsOneGroup{
    my $dir = shift;
    my %low = ();
    my %high = ();
    my $low_file = $dir."/low_confidence_final_results.csv";
    my $high_file = $dir."/high_confidence_final_results.csv";
    my $aux = `ls $dir/*_high.tab`;
    my @files = split("\n", $aux);
    for my $file (@files){
        if($file eq "." or $file eq ".."){
            next;
        }
        open(FILE, $file);
        while(<FILE>){
            chomp($_);
            if($_ =~ /^#/){
                next;
            }
            my @line = split(" ", $_);
            if(defined $high{$line[0]}){
                $high{$line[0]} += 1;
            }
            else{
                $high{$line[0]} = 1;
            }
        }
        close(FILE);
    }
    open(OUT, ">$high_file");
    my %selected_high = ();
    foreach my $key (sort keys %high){
        if($high{$key} > 1){
            $selected_high{$key} = 1;
            print OUT "$key\n";
        }
    }
    close(OUT);
    $aux = `ls $dir/*_low.tab`;
    @files = split("\n", $aux);
    for my $file (@files){
        if($file eq "." or $file eq ".."){
            next;
        }
        open(FILE, $file);
        while(<FILE>){
            chomp($_);
            if($_ =~ /^#/){
                next;
            }
            my @line = split(" ", $_);
            if(defined $low{$line[0]}){
                $low{$line[0]} += 1;
            }
            else{
                $low{$line[0]} = 1;
            }
        }
        close(FILE);
    }
    open(OUT, ">$low_file");
    foreach my $key (sort keys %low){
        if($low{$key} > 1 and !defined $selected_high{$key}){
            print OUT "$key\n";
        }
    }
    close(OUT);

}
