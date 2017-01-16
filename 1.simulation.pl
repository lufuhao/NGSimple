#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw/max min/;
use Cwd;
use FindBin qw($Bin);
use constant USAGE=><<EOH;

SYNOPSIS

perl $0 -i reference.fa --lib_type se|pe|mp --lib_read_length INT --lib_size INT_list --lib_cov INT_list --seed INT --output_dir PATH --output_fmt fa|fq
	
	INT_list: saperated by comma
	PATH: output path/directory
	
VERSION
	20170117

DESCRIPTION
	Simulation reads with specied insert sizes, read length and sequencing coverages;

REQUIREMENT
	PATH: 		mason
	PerlModules: 	GetOpt::Long
			Cwd
			FindBin

OPTIONS
	--help/-h
		Print this help/usage;
	--input|-i <FILE>
		[Msg] Input reference fasta file to simulate the reads;
	--config|-c <File>
		[Opt] Config file specify parameters;
	--lib_type <se|pe|mp>
		[Opt] single-end | paired-end | mate-paired or a mixed 
		list delimited by comma, like: se,pe,pe,mp
	--lib_cov <FLOAT | FLOAT_list>
		[Opt] Sequencing coverage for simulated reads,
		Value: Int or float or list (delimited by comma);
	--lib_read_length <INT|LIST>
		[Opt] Read length of simulated reads;
		Value: Int or list delimited by comma, like: 50 or 20,30,5,1;
	--lib_insertsize <INT|LIST>
		[Opt] Fragmentation loop size for PE and MP;
		Value: Int or list delimited by comma, like: 500 or 500,3000,40000,5000;
	--lib_size_stdev <INT>
		[Opt] Standard deveriation of Fragmentation loop size; Default: 10;
		For 3K, 3000*10%=300 -->  2700<=insert_size<=3000bp;
	--seed <INT>
		[Optional] The seed for Rng in mason, Default: 88;
	--output_dir <PATH>
		[Opt] a NEW folder is recommenred for subsequent anylysis;
		Default: ./;
	--mason_path <PATH/to/mason>
		[Opt] Path to mason if not in PATH;
	--ref_length <INT>
		[Opt] total length of input fasta to calculate number of reads;
		Default: internal calculation by counting characters;
	--verbose
		[Opt] Detailed output for trouble-shooting;
	--version/-v
		[Opt] Print current SCRIPT version;

EXAMPLE
	perl 1.simulation.pl -i my.ref.fasta --lib_type se --lib_cov 10 --seed 88 --output_fmt fq --lib_read_length 100
	perl 1.simulation.pl -i my.ref.fasta --lib_type pe --lib_cov 10 --seed 88 --output_fmt fq --lib_read_length 100 --lib_insertsize 7000
	perl 1.simulation.pl -i my.ref.fasta --lib_type mp --lib_cov 10 --seed 88 --output_fmt fq --lib_read_length 100 --lib_insertsize 7000

AUTHOR
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk
EOH
die USAGE unless @ARGV;



our ($help_message, $input_fasta, $config_file, $seed, $output_dir, $output_format, $fasta_total_length, $masonpath, $verbose, $version);
our ($lib_type, $lib_read_length, $lib_size, $lib_size_stdev, $lib_cov, @lib_types, @lib_read_lengths, @lib_sizes, @lib_covs);

GetOptions (
'help|h!'=> \$help_message,
'reference|i:s' => \$input_fasta,
'config|c:s' => \$config_file,
'lib_insertsize:s' => \$lib_size,
'lib_size_stdev:i' => \$lib_size_stdev,
'lib_cov:i' => \$lib_cov,
'lib_type:s' => \$lib_type,
'lib_read_length:s' => \$lib_read_length,
'seed|s:i'=> \$seed,
'output_dir|o:s'=> \$output_dir,
#'output_fmt:s'=> \$output_format,
'mason_path:s'=> \$masonpath,
'ref_length'=> \$fasta_total_length,
'verbose!'=> \$verbose,
'version|v!'=> \$version,
);
($help_message or $version) and die USAGE;



###Input test##############################
our %RCFcomb=();
our %RCFlib=();
our @lib_names=();
if (defined $config_file and -s $config_file) {
	&ReadConfigureFile($config_file);
}
else {
	if (defined $input_fasta and -s $input_fasta) {
		die "Please input the fasta reference\n";
	}
	else{
		print "Your input fasta file is $input_fasta\n";
	}
	if (defined $lib_type) {
		@lib_types=split(/,/, $lib_type);
		print "Library_type: @lib_types\n";
	}
	else {
		die "Please input the lib_type, or specify that in config_file\n";
	}
	if (defined $lib_read_length) {
		@lib_read_lengths=split(/,/, $lib_read_length);
		print "Library_read_length: @lib_read_lengths\n";
	}
	else {
		die "Please input the read_length or specify that in config_file\n";
	}
	$lib_size=0 unless (defined $lib_size);
	@lib_sizes=split(/,/, $lib_size);
	print "Library_size: @lib_sizes\n";
	if  (defined $lib_cov) {
		@lib_covs=split(/,/, $lib_cov);
		print "Library_coverage: @lib_covs\n";
	}
	else {
		die "Please input the coverage corresponding to the lib_size\n";
	}
}
###test if the arrays get the same lenth
my $max_num_sim=max scalar(@lib_types),scalar(@lib_read_lengths),scalar(@lib_covs);
my $min_num_sim=min scalar(@lib_types),scalar(@lib_read_lengths),scalar(@lib_covs);
die "Empty lib_types, lib_read_lengths or lib_covs\n" if ($min_num_sim<1 or $max_num_sim<1);
&AutoFilArr($max_num_sim, @lib_types);
&AutoFilArr($max_num_sim, @lib_read_lengths);
&AutoFilArr($max_num_sim, @lib_covs);
&AutoFilArr($max_num_sim, @lib_sizes);
unless (scalar(@lib_names)>1) {
	for (my $i=1; $i<=$max_num_sim; $i++) {
		push(@lib_names, "Lib$i");
	}
}
print "Info: A total of $max_num_sim libs are setting to be produced\n" if (defined $verbose);



###test output
$output_dir=getcwd unless (defined $output_dir);
if (-e $output_dir) {
	print "Warning: The output path exists and seems good\n" if (defined $verbose);
}
else {
	mkdir ("$output_dir", 0766) || die "Can not create the output_dir\n";
}
print "The simulated reads will be output to: $output_dir\n" if (defined $verbose);
###correct output format
$output_format="fq";# unless (defined $output_format);
#if (lc($output_format) =~ /(^fasta)|(^fast)|(^fst)|(^fsa)|(^ft)|(^fs)|(^fa)|(^fas)/i) {
#	$output_format='fa';
#}
#elsif (lc($output_format) =~ /(^fastq)|(^fq)/i) {
#	$output_format='fq';
#}
#else {
#	die "Unknown reads format: lc($output_format)\nPlease use one of (fa, fq)\n";
#}
print "The output format was set to: $output_format\n";
###test mason
$masonpath='mason' unless (defined $masonpath);
#`$masonpath illumina -h`;
my $mason_test=&exec_cmd("which $masonpath");
#die "Can not find mason program in \@PATH\nTest_code: $mason_test\nPlease spceify --ref_length PATH/TO/mason if you already installed it\n" if ($mason_test);
###End#######################################



###Default value#############################
$lib_size_stdev=10 unless (defined $lib_size_stdev);
$seed=88 unless (defined $seed);
$fasta_total_length=0 unless (defined $fasta_total_length);
###End#######################################



###Calculate.total.length.of.fasta###########
if ($fasta_total_length == 0) {
	(my $num_seq, $fasta_total_length) = &fastasum($input_fasta);
	print "\n\nThe total length of your fastq file is $fasta_total_length\n\n" if (defined $verbose); ###For test###
}
$fasta_total_length==0 and die "Zero-size reference: $input_fasta\n";
###End#######################################



###Generating the reads######################
my $num_successful_sim=0;
chdir "$output_dir" or die "Can not cd dir $output_dir£º$!\n";
mkdir ("$output_dir/1.fastq", 0766) || die "Can not create directory $output_dir/1.fastq\n";
chdir "$output_dir/1.fastq" or die "Can not cd dir \"$output_dir/1.fastq\"£º$!\n";
for (my $i=0; $i<$max_num_sim; $i++) {
	my @lib=();
	@lib=($input_fasta, $seed, $lib_read_lengths[$i], $output_format, $lib_covs[$i], $lib_types[$i], $lib_sizes[$i], $lib_size_stdev);
#required: 	my ($input, $seed, $read_length, $format, $cov, $type, $insert_size, $insert_size_sd);
	if (&exec_mason(@lib)==0){
		print "The lib with \nlib_types: $lib_types[$i]\ninsert size: $lib_sizes[$i]\nCoverage of $lib_covs[$i], and Readlength: $lib_read_lengths[0]\n finished\n";
		$num_successful_sim++;
	}
}
print "\n\n\nSummary: $num_successful_sim OUT of total $max_num_sim libs successfully finished\n\n\n";
chdir "$output_dir" or die "Can not cd dir $output_dir£º$!\n";
###End##############################



#######################################
###Subfunction#########################
#######################################
###Read parameters from config file
#&ReadConfigureFile($Config_file)
###Need global variables: $input_fasta, @lib_types, @lib_read_lengths, @lib_covs, @lib_sizes, %RCFcomb, @lib_names, %RCFlib
sub ReadConfigureFile {
	my $config_file=shift @_;
	open (CONFIG, "$config_file") || die "Can not open file:$config_file\n";
	my @RCFlibs=();
	%RCFlib=();
	@lib_names=();
	my @RCFcombination=();
	%RCFcomb=();
	my $RCFtest_libs=0;
	my $RCFtest_combination=0;
	my $RCFtest_getfasta=0;
	while (my $line=<CONFIG>) {
		chomp $line;
		next if (($line eq '') or $line=~m/^#|^\*/);
		if ($RCFtest_getfasta==0 and $line=~/^\S+=(\S+)\s*/) {
			$input_fasta=$1;
			$RCFtest_getfasta=1;
		}
		if ($line=~/\@library/i) {
			$RCFtest_libs=1;
			next;
		}
		if ($line=~/\@indice/i) {
			$RCFtest_libs=0;
			$RCFtest_combination=1;
			print "Reading library setting \@library from config_file finished\n";
			next;
		}
		if ($RCFtest_libs==1) {
			@RCFlibs=split(/\t/, $line);
			if (scalar(@RCFlibs)==5) {
				(defined $RCFlibs[0] and $RCFlibs[0] ne '' and $RCFlibs[0] !~ /\s+/) ? push(@lib_names, $RCFlibs[0]) : die "Unknown $RCFlibs[0] for libary $RCFlibs[0] in CONTIG file\n";
				(defined $RCFlibs[1] and $RCFlibs[1] ne '' and $RCFlibs[1]=~/se|pe|mp/i) ?  push(@lib_types, $RCFlibs[1]) : die "Unknown $RCFlibs[1] for library $RCFlibs[0] in CONFIG file\n";
				(defined $RCFlibs[2] and $RCFlibs[2] ne '' and $RCFlibs[2]>0) ? push(@lib_read_lengths, $RCFlibs[2]) : die "Unknown $RCFlibs[2]  for library $RCFlibs[0] in CONFIG file\n";
				(defined $RCFlibs[3] and $RCFlibs[3] ne '' and $RCFlibs[3]>0) ? push(@lib_covs, $RCFlibs[3]) : die "Unknown $RCFlibs[3]  for library $RCFlibs[0] in CONFIG file\n";
				$RCFlibs[4]=0 if ($RCFlibs[1]=~/se/i);
				(defined $RCFlibs[4] and $RCFlibs[3] ne '' and $RCFlibs[3]>=0) ? push(@lib_sizes, $RCFlibs[4]) : die "Unknown $RCFlibs[4]  for library $RCFlibs[0] in CONFIG file\n";
###Format: @{$RCFlib{libname}}=(lib_type, readlength, coverage, insersize);
				@{$RCFlib{$RCFlibs[0]}}=($RCFlibs[1], $RCFlibs[2], $RCFlibs[3], $RCFlibs[4]);
				print "Library $RCFlibs[0] Accepted: $RCFlibs[1], $RCFlibs[2], $RCFlibs[3], $RCFlibs[4]\n";
			}
			else {
				next;
			}
		}
		if ($RCFtest_combination==1) {
			@RCFcombination=split(/\t/, $line);
			if (scalar(@RCFcombination)>=2 and scalar(@RCFcombination)<=11) {
				for (my $j=0; $j<scalar(@RCFcombination)-1; $j++) {
					die "Unknown conbination $RCFcombination[0]\n" unless (exists $RCFlib{$RCFcombination[$j+1]});
					@{${$RCFcomb{$RCFcombination[0]}}[$j]}=@{$RCFlib{$RCFcombination[$j+1]}};
				}
			}
		}
	}
	close CONFIG;
	my @RCFindices=keys %RCFcomb;
	foreach my $RCFindex (@RCFindices) {
		my $RCF_num_libs=scalar(@{$RCFcomb{$RCFindex}});
		print "Index $RCFindex Accepted: \n";
		for (my $i=0; $i<$RCF_num_libs; $i++) {
			print "\t@{${$RCFcomb{$RCFindex}}[$i]}\n";
		}
	}
	print "Reading library combination \@indices from Config_file finished\n";
}



##Add full path to a file name
#&AddFullPath($file, 0/1);
sub AddFullPath {
	my ($file, $test_exist) = @_;
	if (ref($file) eq "ARRAY"){
		for (my $i=0;$i<scalar(@$file);$i++){
			my $filename = $file->[$i];
			if ($test_exist && (! -e $filename)) {
				die "Error, cannot locate file: $filename";
			}
			$file->[$i] = &create_full_path($filename);
		}
		return @$file;
	}
	else {
		if ($test_exist && ! -e $file) {
			die "Error, cannot locate file: $file";
		}
		use Cwd;
		my $cwd=getcwd;
		no Cwd;
		if ($file !~ m/^\//) {
			$file=$cwd.'/'.$file;
		}
		return($file);
	}
}



sub AutoFilArr {
	my ($num, @array)=@_;
	my $constant='';
	$constant=$array[-1];
	@array=map {defined($_)?$_:0} @array; ###reset undefined elements to 0
	for (my $i=scalar(@array); $i<$num; $i++) {
		push @array, $constant;
	}
	return @array;
}



sub fastasum {
	my $FSfile=shift @_;
	open (INPUT, "$FSfile") || die "Can not open fasta files specified: $FSfile\n";
	my $FSnum_seq=0;
	my $FStotal_length=0;
	my $FSseqid='';
	my $FSline_num=0;
	while (my $FSline=<INPUT>) {
		chomp $FSline;
		my $FSline_num++;
		if ($FSline =~/^>/) {
			$FSnum_seq++;
			$FSseqid=$FSline;
		}
		else {
			$FStotal_length+=length($FSline);
			if ($FSline=~/[^atcgATCGnN]/) {
				print "Warning: Non-ATCGatcgNn letters at line $FSline_num in sequence $FSseqid in file $FSfile\n";
			}
		}
	}
	close INPUT;
	return ($FSnum_seq, $FStotal_length);
}



###Design Mason input###############################################
###Mason simulate reads###################################################
###$mason illumina -s 88 -N 200000 -o is3kb_rl100.fastq -sq -i -ll 3000 -le 300 -lu -mp -n 100 -rn 2 reference.fasta
###mason illumina --help
###      -s	--seed INT	The seed for Rng. Default: 0.
###      -N	--num-reads INT	Number of reads (or mate pairs) to simulate. Default: 1000. 
###      -o	--output-file FILE	Write results to PARAM.fasta file instead of SEQUENCE.reads.fasta. Valid filetypes are fasta,fa,fastq,and,fq
###-#sq--#simulate-#qualities#Simulate qualities, generate #FASTQ instead of #FASTA;
###      -i	--include-read-information	Include additional read information in reads file.
###      -ll	--library-length-mean NUM	Mate-pair mean library length. Default: 1000. 
###      -le	--library-length-error NUM	Mate-pair library tolerance. Default: 100. 
###      -lu	--library-length-uniform	Mate-pair library length is uniform. 
###      -mp	--mate-pairs	Enable mate pair simulation. 
###      -n	--read-length NUM	The length of the reads to simulate. All resulting reads will have the same length. Default: 36.
###      -rn	--read-naming NUM	Read naming scheme in FASTQ/FASTA files. See section on read naming below. In range [0..2]. Default: 0. 
###      		READ NAMING SCHEME 
###      		Note that the suffixes will not appear in the SAM file. In the SAM format, 
###      		the FLAG field is used to mark the leftmost/rightmost fragment/read. 
###      		    0     No suffix, left-end and right-end read will be named 
###      		          'example.fastq.0000', for example.
###      		    1     Add zero-based slash-suffix, i.e. names will end on '/0' and '/1' 
###     		    2     Add one-based slash-suffix, i.e. names will end on '/1' and '/2' 
###      -rnp	--read-name-prefix STR	Read name prefix. Default is output file name.
####################################################################################
#output: $output_dir/1.fastq/Sim.$type.$insert_size.$cov.$read_length.$mason_outputformat
sub exec_mason {
	my ($input, $seed, $read_length, $format, $cov, $type, $insert_size, $insert_size_sd)=@_;
	my $cmd='';
	my $num_reads=0;
	my $insert_error=int($insert_size*$insert_size_sd/100);
	my $mason_cmd=$masonpath." illumina -s $seed -n $read_length ";
#decide read format
	my $mason_outputformat='';
	if ($format eq 'fa') {
		$mason_outputformat='fa';
	}
	elsif ($format eq 'fq') {
		$mason_outputformat='fq';
		$mason_cmd.=' -sq ';
	}
#decide libtype
	if ($type eq 'se') {
		$num_reads=int($fasta_total_length*$cov/$read_length);
	}
	elsif ($type eq 'pe' or $type eq 'mp') {
		$num_reads=int($fasta_total_length*$cov/(2*$read_length));
		$mason_cmd.=" --library-length-mean $insert_size --library-length-error $insert_error --read-naming 2 -rnp Sim$insert_size.$cov --mate-pairs --library-length-uniform ";
	}
#Read number control due to Mason will refuse to work when read number is huge
	if ($num_reads<=30000000) {
		my $EMoutput='';
		$EMoutput="./Sim.$type.$insert_size.$cov.$read_length.$mason_outputformat";
		unlink($EMoutput) if (-e $EMoutput);
		$mason_cmd.=" -o $EMoutput -N $num_reads $input ";
		&exec_cmd($mason_cmd);
	}
	else {
		my $total_reads=$num_reads;
		my $indv_numread=30000000;
		my $num_split=int($num_reads/30000000)+1;
		for (my $i=1; $i<=$num_split;$i++) {
			my $lane='L'.$i;
			my $EMoutput='';
			$EMoutput="./Sim.$type.$insert_size.$cov.$read_length.$lane.$mason_outputformat";
			unlink ($EMoutput) if (-e $EMoutput);
			if ($total_reads>$indv_numread) {
				$mason_cmd.=" -o $EMoutput -N 30000000 $input";
			}
			else {
				$mason_cmd.=" -o $EMoutput -N $total_reads $input";
			}
			&exec_cmd($mason_cmd);
			$total_reads-=$indv_numread;
		}
		$cmd="cat ./Sim.$type.$insert_size.$cov.$read_length.L*.$mason_outputformat > ./Sim.$type.$insert_size.$cov.$read_length.$mason_outputformat";
		&exec_cmd($cmd);
		unlink("./Sim.$type.$insert_size.$cov.$read_length.L*.$mason_outputformat");
	}
	return 0;
}



sub mytime {
	my($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime();
	$year += 1900;
	$mon  += 1;
	my $btime = sprintf("%04d%02d%02d %02d:%02d:%02d",$year,$mon,$day,$hour,$min,$sec);
	return $btime;
}



sub exec_cmd {
	my ($cmd)=@_;
	print &mytime()."CMD: $cmd\n";
	my $start_time=time();
	my $return_code=system($cmd);
	my $end_time=time();
#	if ($return_code == -1) {
#		print ¡°failed to execute: $!\n¡±;
#	}
#	elsif ($return_code & 127) {
#		printf ¡°child died with signal %d, %s coredump\n¡±, ($return_code & 127),  ($return_code & 128) ? ¡®with¡¯ : ¡®without¡¯;
#	}
#	else {
#		printf ¡°child exited with value %d\n¡±, $return_code >> 8;
#	}
	if ($return_code) {
		die "Error, cmd: $cmd died with ReturnCode $return_code\n";
		return $return_code;
	}
	else {
		print "Finished command: $cmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n";
		return $return_code;
	}
}
