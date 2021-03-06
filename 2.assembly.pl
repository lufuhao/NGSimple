#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw/max min/;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
use constant USAGE=><<EOH;

SYNOPSIS

perl $0 --config CONFIG --output_dir ./
	
VERSION
	LUFUHAO20170206

DESCRIPTION
	Assemble the simulated reads under output_dir/1.Fastq by 
	setting in config file, choose the Kmer with largest n50, 
	under under output_dir/2.Assembly, evaluate the asseblies,
	And visualize the assembly against reference under 
	output_dir/3.summary

REQUIREMENT
	PATH: 		java, velvetg, velveth, nucmer, 
	PerlModules: 	Getopt::Long
			List::Util
			Cwd
			FindBin

OPTIONS
	--help|-h
		Print this help/usage;
	--config|-c <File>
		[Msg] Configure file;
	--reference|-r <File>
		[Opt] Reference fasta file or specify in config;
	--threads|-t <INT>
		[Opt] Number of running threads;
	--output_dir|-d <PATH>
		[Opt] Path containing simulated reads: '1.Fastq';
		Default: ./
	--ploteachK
		[Opt] Whether run mummplot for each K
	--path_java </Path/to/java>
		[Opt] path to java if Not in PATH;
	--path_trimmomatic </Path/to/trimmomatic.jar>
		[Opt] path to trimmomatic-verion.jar; Defaults: 
		PROGRAM_ROOT/utils/Trimmomatic/v0.32/x86_64/trimmomatic-0.32.jar
	--minlength|-l <INT>
		[Opt] Minimum read length when trimming to keep a read;
		Default: 36;
	--path_velvetg </Path/to/velvetg>
		[Opt] Path to velvetg if NOT in PATH;
	--path_velveth </PATH/to/velveth>
		[Opt] Path to velveth if NOT in velveth;
	--path_nucmer </Path/to/nucmer>
		[Opt] Path to nucmer in Mummer if NOT in PATH;
	--path_deltafilter <PATH/to/delta-filter>
		[Opt] Path to delta-filter if NOT in PATH;
	--path_mummerplot </Path/to/mummerplot>
		[Opt] Path to mummerplot if NOT in PATH;
	--verbose
		[Opt] Detailed output for trouble-shooting;
	--version|-v
		[Opt] Print current SCRIPT version;

EXAMPLE
	perl --config ./CONFIG --output_dir ./

AUTHOR
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk
EOH
die USAGE unless @ARGV;



my ($help_message, $verbose, $version);
my ($config_file, $output_dir, $test_run_stepwiseK, $num_threads);
my ($path_java, $path_trimmomatic, $trimmo_minlen);
my ($input_fasta, $path_velvetg, $path_velveth);
my ($path_nucmer, $path_deltafilter, $path_mummerplot);
my $quality_score_format;

my $mummer_each_K=0;### 1: run mummplot for each K; 0 just for best K
GetOptions (
'help|h!'=> \$help_message,
'config|c:s' => \$config_file,
'reference|r:s' => \$input_fasta,
'threads|t:i' => \$num_threads,
'output_dir|d:s' => \$output_dir,
'path_java:s' => \$path_java,
'path_trimmomatic:s' => \$path_trimmomatic,
'minlength|l:i' => \$trimmo_minlen,
'ploteachK!' => \$mummer_each_K,
'path_velvetg:s' => \$path_velvetg,
'path_velveth:s' => \$path_velveth,
'path_nucmer:s' => \$path_nucmer,
'path_deltafilter:s' => \$path_deltafilter,
'path_mummerplot:s' => \$path_mummerplot,
'verbose!'=> \$verbose,
'version|v!'=> \$version,
);
($help_message or $version) and die USAGE;



###Defaults and test###################
my (@lib_names, @lib_types, @lib_read_lengths, @lib_covs, @lib_sizes, %RCFlib, %RCFcomb, $MaxVelvetCategories);
$num_threads=1 unless (defined $num_threads);
my $experence_factor=0.3;
my $cmd='';
$path_java='java' unless (defined $path_java);
$path_trimmomatic="$Bin/utils/Trimmomatic/v0.32/x86_64/trimmomatic-0.32.jar" unless (defined $path_trimmomatic);
$quality_score_format='phred33' unless (defined $quality_score_format);
die "Error: Quality format accept \'phred33\' or \'phred64\' only\n" unless ($quality_score_format eq 'phred33' or $quality_score_format eq 'phred64');
$test_run_stepwiseK=1 unless (defined $test_run_stepwiseK);
$trimmo_minlen=36;
$MaxVelvetCategories=10;
$path_velveth='velveth' unless (defined $path_velveth);
$path_velvetg='velvetg' unless (defined $path_velvetg);
$path_nucmer='nucmer' unless (defined $path_nucmer);
$path_deltafilter='delta-filter' unless (defined $path_deltafilter);
$path_mummerplot='mummerplot' unless (defined $path_mummerplot);
my $test_deletetempKfolder=1;
my @longestcontigfiles=();
my @allcontigfile=();

###test output
my $run_dir=getcwd;
$output_dir=$run_dir unless (defined $output_dir);
if (-e $output_dir) {
	print "Warning: The output path exists and seems good\n" if (defined $verbose);
}
else {
	mkdir ("$output_dir", 0766) || die "Can not create the output_dir\n";
}
print "The simulated reads will be output to: $output_dir\n" if (defined $verbose);

###End#################################



###Main############################
if (-d "$output_dir/2.assembly") {
	&DeletePath("$output_dir/2.assembly");
}
mkdir ("$output_dir/2.assembly", 0766) || die "Error: can not create directory:$output_dir/2.assembly\n";
&ReadConfigureFile($config_file);
die "Error: Please specify your referece sequence in CONFIG file\n" unless (defined $input_fasta and $input_fasta ne '');
foreach my $AssemblyIndex (sort keys %RCFcomb) {
	print "Assembling Index: $AssemblyIndex\n";
	unless (&SettingVelvet($AssemblyIndex)) {
		print STDERR "Warnings: Failed assembly: $AssemblyIndex\n";
		next;
	}
}
unless (&AssemblySummary()) {
	die "Error: collect summary failed, You need manually check\n";
}
###End#################################



#######################################
###Subfunction#########################
#######################################
###Step1: Read parameters from config file
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
			chomp $input_fasta;
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
			if (scalar(@RCFcombination)>=2 and scalar(@RCFcombination)<=($MaxVelvetCategories+1)) {
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
	
	return 1;
}



###Step2: Set velvet
#&SetVelvet($index)
#global variables: %RCFcomb, $test_run_stepwiseK, $input_fasta, $experence_factor, @longestcontigfiles=(); @allcontigfile=();
sub SettingVelvet {
	my $SVindex=shift @_;
	chdir "$output_dir" or die "Velvet error: can not cd dir $output_dir：$!\n";
	my $newdir="$output_dir/2.assembly/$SVindex";
	&exec_cmd("rm -rf $newdir") if (-d $newdir);
	mkdir ($newdir, 0766) || die "Velvet error: can not create directory $newdir\n";
	chdir "$newdir" or die "Vetvet error: can not cd dir $newdir：$!\n";
	my @SVseq_files=();
	my $SVreadset='';
	my $SVvelvetg_set='';
	my $SV_num_lib=1;
	my $SV_num_lib2=1;
	my $SVtotal_coverage=0;
	my $SVmin_readlength=10000;
	for (my $SVi=0; $SVi<scalar(@{$RCFcomb{$SVindex}}); $SVi++) {
		if (${${$RCFcomb{$SVindex}}[$SVi]}[0] eq 'se') {
			my $SVfastq_file='';
			($SVfastq_file)=&ReadFastqFiles(${${$RCFcomb{$SVindex}}[$SVi]}[0], ${${$RCFcomb{$SVindex}}[$SVi]}[3], ${${$RCFcomb{$SVindex}}[$SVi]}[2], ${${$RCFcomb{$SVindex}}[$SVi]}[1]);
			push (@SVseq_files, $SVfastq_file);
			$SVreadset.=" -short -fastq $SVfastq_file";
			$SV_num_lib2++;
		}
		elsif (${${$RCFcomb{$SVindex}}[$SVi]}[0] eq 'pe' or ${${$RCFcomb{$SVindex}}[$SVi]}[0] eq 'mp') {
			my ($SVfastq_R1, $SVfastq_R2)=('', '');
			($SVfastq_R1, $SVfastq_R2)=&ReadFastqFiles(${${$RCFcomb{$SVindex}}[$SVi]}[0], ${${$RCFcomb{$SVindex}}[$SVi]}[3], ${${$RCFcomb{$SVindex}}[$SVi]}[2], ${${$RCFcomb{$SVindex}}[$SVi]}[1]);
			push (@SVseq_files, $SVfastq_R1);
			push (@SVseq_files, $SVfastq_R2);
			if ($SV_num_lib==1) {
				$SVreadset.=" -shortPaired -fastq -separate $SVfastq_R1 $SVfastq_R2 ";
				$SVvelvetg_set.=" -ins_length ${${$RCFcomb{$SVindex}}[$SVi]}[3] -ins_length_sd ". ${${$RCFcomb{$SVindex}}[$SVi]}[3]/5 . ' ';
			}
			else {
				$SVreadset.=" -shortPaired$SV_num_lib -fastq -separate $SVfastq_R1 $SVfastq_R2 ";
				$SVvelvetg_set.=" -ins_length$SV_num_lib2 ${${$RCFcomb{$SVindex}}[$SVi]}[3] -ins_length$SV_num_lib2".'_sd '. int(${${$RCFcomb{$SVindex}}[$SVi]}[3]/5) . ' ';
			}
			$SV_num_lib++; $SV_num_lib2++;
		}
		else {
			die "Unknown libary type: ${${$RCFcomb{$SVindex}}[$SVi]}[0]\n";
		}
		$SVmin_readlength=min $SVmin_readlength,${${$RCFcomb{$SVindex}}[$SVi]}[1];
		$SVtotal_coverage+=${${$RCFcomb{$SVindex}}[$SVi]}[2];
	}
	my $SVgenome_size=&FastaSum($input_fasta);
	unless ($SVgenome_size) {
		print STDERR "Error: targetSequence length 0\n";
		return 0;
	}
	my $SVtargetKmer=$SVtotal_coverage*$experence_factor;
	print "\n\nReference size $SVgenome_size\nExpect K-mer:$SVtargetKmer\nAssembly fies: \n";
	map {print $_."\n"} @SVseq_files;
	my $SVguess_beskK=&GuessBestKmer($SVgenome_size, $SVtargetKmer, @SVseq_files);
	if ($SVguess_beskK >0) {
		print "Info: Guess bast K-mer: $SVguess_beskK\n\n";
	}
	else {
		print STDERR "Error: failed to guess bestK\n";
		return 0;
	}
	
	unless (chdir "$newdir") {
		print STDERR "Error: Vetvet error: can not cd dir $newdir：$!\n";
		return 0;
	}
	if ($test_run_stepwiseK==0) {
		unless (&RunVelvet($SVguess_beskK, $SVreadset, $SVvelvetg_set, $SVindex)) {
			print STDERR "Error: Velvet running failed\n";
			return 0;
		}
		link("$newdir/K$SVguess_beskK/contigs.fa", "$SVindex.K$SVguess_beskK.contigs.fa");
		link("$newdir/K$SVguess_beskK/longest_contig.fa", "$SVindex.K$SVguess_beskK.LongestContig.fa");
		if (-s "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.AllContigs.png") {
			link("$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.AllContigs.png", "$SVindex.K$SVguess_beskK.AllContigs.png");
		}
		if (-s "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.AllContigs.ps") {
			link("$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.AllContigs.ps", "$SVindex.K$SVguess_beskK.AllContigs.ps");
		}
		if (-s "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.LongestContig.png") {
			link("$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.LongestContig.png", "$SVindex.K$SVguess_beskK.LongestContig.png");
		}
		if (-s "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.LongestContig.ps") {
			link("$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.LongestContig.ps", "$SVindex.K$SVguess_beskK.LongestContig.ps");
		}
		
		push (@allcontigfile, "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.AllContigs.fplot") if (-s "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.AllContigs.fplot");
		push (@allcontigfile, "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.AllContigs.rplot") if (-s "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.AllContigs.rplot");
		push (@allcontigfile, "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.AllContigs.gp") if (-s "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.AllContigs.gp");
		
		
		push (@longestcontigfiles, "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.LongestContig.fplot") if (-s "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.LongestContig.fplot");
		push (@longestcontigfiles, "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.LongestContig.rplot") if (-s "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.LongestContig.rplot");
		push (@longestcontigfiles, "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.LongestContig.gp") if (-s "$newdir/K$SVguess_beskK/$SVindex.K$SVguess_beskK.LongestContig.gp");
		
	}
	else {
		my $SVrange=int($SVguess_beskK/5);
		$SVrange+=1 if ($SVrange%2==1);
		my $SVbestKmer=&SelectBestAssembly($SVguess_beskK-$SVrange, $SVguess_beskK+$SVrange, 8, $SVreadset, $SVvelvetg_set, $SVindex);
		unless ($SVbestKmer>0) {
			print STDERR "Error: Failed to select best K ", ($SVguess_beskK-$SVrange), " to ", ($SVguess_beskK+$SVrange), "\n";
			return 0;
		}
		my $SVbestKmer2=&SelectBestAssembly($SVbestKmer-4, $SVbestKmer+4, 4, $SVreadset, $SVvelvetg_set, $SVindex);
		unless ($SVbestKmer2>0) {
			print STDERR "Error: Failed to select best K ", ($SVguess_beskK-4), " to ", ($SVguess_beskK+4), "\n";
			return 0;
		}
		my $SVbestKmer3=&SelectBestAssembly($SVbestKmer-2, $SVbestKmer+2, 2, $SVreadset, $SVvelvetg_set, $SVindex);
		unless ($SVbestKmer3>0) {
			print STDERR "Error: Failed to select best K ", ($SVguess_beskK-2), " to ", ($SVguess_beskK+2), "\n";
			return 0;
		}
		if (-s "$newdir/K$SVbestKmer3/contigs.fa" and -s "$newdir/K$SVbestKmer3/longest_contig.fa") {
			link("$newdir/K$SVbestKmer3/contigs.fa", "$SVindex.K$SVbestKmer3.contigs.fa");
			link("$newdir/K$SVbestKmer3/longest_contig.fa", "$SVindex.K$SVbestKmer3.LongestContig.fa");
			unless ($mummer_each_K) {
				unless (&MumStats($input_fasta, "$newdir/K$SVbestKmer3/contigs.fa", "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.AllContigs")) {
					print STDERR "Error: MUMmerplot running failed: AllContigs $SVindex K$SVbestKmer3\n";
					return 0;
				}
				unless (&MumStats($input_fasta, "$newdir/K$SVbestKmer3/longest_contig.fa", "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.LongestContig")) {
					print STDERR "Error: MUMmerplot running failed: LongestContig $SVindex K$SVbestKmer3\n";
					return 0;
				}
			}
		}
		else {
			print STDERR "Error: Failed to find assemblies: $newdir/K$SVbestKmer3/contigs.fa and $newdir/K$SVbestKmer3/longest_contig.fa\n";
			return 0;
		}
		if (-s "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.AllContigs.png") {
			link("$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.AllContigs.png", "$SVindex.K$SVbestKmer3.AllContigs.png");
		}
		if (-s "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.AllContigs.ps") {
			link("$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.AllContigs.ps", "$SVindex.K$SVbestKmer3.AllContigs.ps");
		}
		if (-s "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.LongestContig.png") {
			link("$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.LongestContig.png", "$SVindex.K$SVbestKmer3.LongestContig.png");
		}
		if (-s "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.LongestContig.ps") {
			link("$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.LongestContig.ps", "$SVindex.K$SVbestKmer3.LongestContig.ps");
		}
		
		push (@allcontigfile, "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.AllContigs.fplot") if (-s "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.AllContigs.fplot");
		push (@allcontigfile, "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.AllContigs.rplot") if (-s "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.AllContigs.rplot");
		push (@allcontigfile, "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.AllContigs.gp") if (-s "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.AllContigs.gp");
		
		
		push (@longestcontigfiles, "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.LongestContig.fplot") if (-s "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.LongestContig.fplot");
		push (@longestcontigfiles, "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.LongestContig.rplot") if (-s "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.LongestContig.rplot");
		push (@longestcontigfiles, "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.LongestContig.gp") if (-s "$newdir/K$SVbestKmer3/$SVindex.K$SVbestKmer3.LongestContig.gp");
		
		if ($test_deletetempKfolder) {
			my @SVfolders=<$newdir/*>;
			foreach my $SVindfolder (@SVfolders) {
				if (-d $SVindfolder and ($SVindfolder ne "$newdir/K$SVbestKmer3")) {
					unless (&exec_cmd("rm -rf $SVindfolder > /dev/null 2>&1")) {
						print STDERR "Warnings: Unable to delete path: $SVindfolder";
					}
				}
			}
		}
	}
	
	unless (chdir "$output_dir") {
		print STDERR "Warnings: Velvet error: can not cd dir $output_dir：$!\n";
	}
	
	return 1;
}



###ReadFastqFiles according to the setting
#&ReadFastqFiles($lib_type, $lib_insertsize, $coverage, $read_length)
#global variables: $output_dir
sub ReadFastqFiles {
	my ($RFFtype, $RFFis, $RFFcov, $RFFrl)=@_;
	my $RFFfile_basename="$output_dir/1.fastq/Sim.$RFFtype.$RFFis.$RFFcov.$RFFrl";
	my $RFFnum_files=0;
	if ($RFFtype eq 'se') {
		unless (-s $RFFfile_basename.'_trim.fq') {
			&Trimmomatic($RFFfile_basename.'.fq');
		}
		if (-s $RFFfile_basename.'_trim.fq'){
			return ($RFFfile_basename.'_trim.fq', 1);
		}
	}
	elsif (($RFFtype eq 'pe') or ($RFFtype eq 'mp')) {
		unless (-s $RFFfile_basename.'_1_trim.fq' and -s $RFFfile_basename.'_2_trim.fq') {
			unlink ($RFFfile_basename.'_1_trim.fq') if (-s $RFFfile_basename.'_1_trim.fq');
			unlink ($RFFfile_basename.'_2_trim.fq') if (-s $RFFfile_basename.'_2_trim.fq');
			if (-s $RFFfile_basename.'_1.fq' and -s $RFFfile_basename.'_2.fq') {
				&Trimmomatic($RFFfile_basename.'_1.fq', $RFFfile_basename.'_2.fq', 2);
			}else {
				die "Error: Can not find Fq files: ".$RFFfile_basename."_1.fq and ".$RFFfile_basename."_2.fq\n";
			}
			
		}
		if (-s $RFFfile_basename.'_1_trim.fq' and -s $RFFfile_basename.'_2_trim.fq') {
			return ($RFFfile_basename.'_1_trim.fq', $RFFfile_basename.'_2_trim.fq')
		}
	}
#	my @RFFfiles=glob "$RFFfile_basename*fq";
#	if ($RFFtype eq 'se') {
#		die "Read fastq error: Can not detect fastq files correctly\n" unless (scalar(@RFFfiles)==1);
#	}
#	elsif (($RFFtype eq 'pe') or ($RFFtype eq 'mp')) {
#		die "Read fastq error: Can not detect fastq files correctly\n" unless (scalar(@RFFfiles)==2);
#	}
#	else {
#		die "Read fastq error: Please do check the $output_dir/1.fastq\n";
#	}
#	foreach my $RFFfastq (@RFFfiles) {
#		die "Detected fastq file is empty: $RFFfastq\n" unless (-s $RFFfastq);
#	}
#	return @RFFfiles;
}



###Trimmomatic trim and clean the reads
#&Trimmomatic(fastqR1, R2, $num_files)
#global varables: $output_dir, $path_java, $path_trimmomatic, $num_threads, $quality_score_format, $trimmo_minlen
sub Trimmomatic {
	my ($TMfastq1, $TMfastq2, $TMnum_files)=@_;
	print "\n\n\n###Trimmomatic BEGIN###\n";
	chdir "$output_dir" or die "Trimmomatic error: can not cd dir $output_dir：$!\n";
	my $newdir=$output_dir.'/1.fastq';
	chdir "$newdir" or die "Trimmomatic error: can not cd dir $newdir：$!\n";
	$cmd='';
	die "Trimmomatic error: Can not find file: $TMfastq1\n" unless (-s $TMfastq1);
	if ($TMnum_files==2) {
		die "Trimmomatic error: Can not find file: $TMfastq2\n" unless (-s $TMfastq2);
	}
	my $TMfile_setting='';
	my $TMfq1_output='';
	my $TMfq2_output='';
	if ($TMnum_files==2) {
		my $TMfq1_basename=&RetrvNoExt($TMfastq1); $TMfq1_output="$output_dir/1.fastq/$TMfq1_basename"."_trim.fq";
		my $TMfq2_basename=&RetrvNoExt($TMfastq2); $TMfq2_output="$output_dir/1.fastq/$TMfq2_basename"."_trim.fq";
		if (-s $TMfq1_output and -s $TMfq2_output) {
			return ($TMfq1_output, $TMfq2_output);
		}
		else {
			unlink($TMfq1_output) if (-s $TMfq1_output);
			unlink($TMfq2_output) if (-s $TMfq2_output);
		}
		$TMfile_setting=" $TMfastq1 $TMfastq2 $TMfq1_output $TMfq1_basename".'_unpaired.fq'." $TMfq2_output $TMfq2_basename".'unpaired.fq ';
		$cmd="$path_java -jar $path_trimmomatic PE";
	} elsif ($TMnum_files==1) {
		my $TMfastq_basename=&RetrvNoExt($TMfastq1); 
		$TMfq1_output=$output_dir.'/1.fastq/'.$TMfastq_basename.'_trim.fq';
		if (-s $TMfq1_output) {
			return ($TMfq1_output);
		}
		$TMfile_setting=" $TMfastq1 $TMfq1_output ";
		$cmd="$path_java -jar $path_trimmomatic SE"
	}
	$cmd.=" -threads $num_threads ";
	if ($quality_score_format eq 'phred33') {
		$cmd.=' -phred33 ';
	}
	elsif ($quality_score_format eq 'phred64') {
		$cmd.=' -phred64 ';
	}
	else {
		die "Trimmomatic error: Uknown quality score format $quality_score_format\n";
	}
	$cmd.=$TMfile_setting."HEADCROP:8 LEADING:18 TRAILING:18 SLIDINGWINDOW:4:18 MINLEN:$trimmo_minlen";
	unless (&exec_cmd($cmd)) {
		die "Error: trimmatic failed\n";
	}
	chdir "$output_dir" or die "Trimmomatic error: can not cd dir $output_dir：$!\n";
	print "###Trimmomatic END###\n\n\n";
	if ($TMfq2_output ne '' and -s $TMfq2_output) {
		return ($TMfq1_output, $TMfq2_output);
	} elsif ($TMfq2_output eq '') {
		return ($TMfq1_output);
	}
}



###GuessKmer for fastq file
#&GuessBestKmer($Genome_size, $Target_kmer, @fastq_files)
sub GuessBestKmer {
	print "Guessing Kmer...\n";
	my ($GBKgenome_size, $GBK_targetKmer, @GBKseq_files)=@_;
	my %kmer_count=();
	foreach my $GBKseq_ind_file (@GBKseq_files) {
		open (INPUT, $GBKseq_ind_file) || die "Can not find file: $GBKseq_ind_file\n";
		while (my $line=<INPUT>) {
			chomp $line;
			if ($line=~m/^>/) {
				$line=<INPUT>;
				chomp $line;
				$kmer_count{length($line)}++;
			}
			if ($line=~m/^@/) {
				$line=<INPUT>;
				chomp $line;
				$kmer_count{length($line)}++;
				<INPUT>;<INPUT>;
			}
		}
		close INPUT;
	}
	my @GBKread_length=sort {$a<=>$b} keys %kmer_count;
	if (scalar(@GBKread_length)<1) {
		print STDERR "Error: Can not find any read\n";
		return 0;
	}
	my $GBKbest_kmer_cov=0;
	my $GBKbestK=0;
	print "K\tNum_kmers\tKmer_cov\tExp_Kmer_cov\n" if (defined $verbose);
	for (my $i=1; $i<=$GBKread_length[-1]; $i+=2) {
		my $num_kmers=0;
		foreach my $GBKkmer (@GBKread_length) {
			next if ($GBKkmer<$i);
			$num_kmers+=($GBKkmer-$i+1)*$kmer_count{$GBKkmer};
		}
		my $GBKKmer_cov=$num_kmers/$GBKgenome_size;
		if (abs($GBKKmer_cov-$GBK_targetKmer)<abs($GBKbest_kmer_cov-$GBK_targetKmer)) {
			$GBKbest_kmer_cov=$GBKKmer_cov;
			$GBKbestK=$i;
		}
		print "$i\t$num_kmers\t$GBKKmer_cov\t$GBK_targetKmer\n" if (defined $verbose);
	}
	($GBKbestK>1) ? (return $GBKbestK) : (return 0);
}



###select largest n50 for ./K*/contigs.fa
#&SelectBestAssembly(min, max, step, $SVreadset, $SVvelvetg_set, $SBAindex)
sub SelectBestAssembly {
	my ($SBAminK, $SBAmaxK, $SBAstep, $SBAreadset, $SBAvelvetg_set, $SBAindex)=@_;
	for (my $SBAj=$SBAminK; $SBAj<=$SBAmaxK; $SBAj+=$SBAstep) {
		unless (-d "./K$SBAj") {
			unless (&RunVelvet($SBAj, $SBAreadset, $SBAvelvetg_set, $SBAindex)) {
				print STDERR "Error: Velvet runing failed\n";
				return 0;
			}
		}
	}
	my $SBAcur_dir=getcwd;
	my $SBAbestK=0;
	my $SBAmax_n50=0;
	my @SBAfiles=glob "$SBAcur_dir/K*/contigs.fa";
	unless (scalar(@SBAfiles)>0) {
		print STDERR "Error: Can not find any $SBAcur_dir/K*/contigs.fa velvet assemblies\n";
		return 0;
	}
	my @SBAKs=();
	foreach my $SBAcontig_file (@SBAfiles) {
		my $SBAind_n50=&ContigsStats($SBAcontig_file, 0, 'n50');
		if ($SBAcontig_file=~m/\/K(\d{1,3})\//) {
			my $SBAcurk=0;
			$SBAcurk=$1;
			print "K$SBAcurk n50: $SBAind_n50\n";
			print "K$SBAcurk stats: ". join("\t", &ContigsStats($SBAcontig_file, 0, 'stats')) ."\n";
			push (@SBAKs, $SBAcurk);
			if ($SBAind_n50>$SBAmax_n50) {
				$SBAmax_n50=$SBAind_n50;
				$SBAbestK=$SBAcurk;
			}
			print "BestK selected: $SBAbestK\n";
		}
	}
	@SBAKs=sort {$a<=>$b} @SBAKs;
	if (defined $SBAbestK and $SBAbestK >0) {
		if ($SBAbestK==$SBAKs[0]) {
#			&SelectBestAssembly($SBAKs[0]-2*$SBAstep, $SBAKs[0]+$SBAstep, $SBAstep, $SBAreadset, $SBAvelvetg_set, $SBAindex);
			&SelectBestAssembly($SBAKs[0]-$SBAstep, $SBAKs[0]+$SBAstep, $SBAstep, $SBAreadset, $SBAvelvetg_set, $SBAindex);
		}
		elsif ($SBAbestK==$SBAKs[-1]) {
#			&SelectBestAssembly($SBAKs[-1]-$SBAstep, $SBAKs[-1]+2*$SBAstep, $SBAstep, $SBAreadset, $SBAvelvetg_set, $SBAindex);
			&SelectBestAssembly($SBAKs[-1]-$SBAstep, $SBAKs[-1]+$SBAstep, $SBAstep, $SBAreadset, $SBAvelvetg_set, $SBAindex);
		}
		elsif ($SBAbestK>$SBAKs[0] and $SBAbestK<$SBAKs[-1]) {
			return $SBAbestK;
		}
		else {
			print STDERR "Error: SelectBestAssembly error\n";
			return 0;
		}
	}
	else {
		print STDERR "Error: Can not calcaulate best Kmer between $SBAminK - $SBAmaxK with a step $SBAstep dor read set \'$SBAreadset\'\n";
		return 0;
	}
}



###RunVelvet
#&RunVelvet($k, $readset, $velvetg_set)
#global variables: $path_velveth, $path_velvetg
sub RunVelvet {
	my ($RVk, $RVreadset, $RVvelvetg_set, $RVindex)=@_;
	my $velveth_cmd='';
	my $velvetg_cmd='';
	my $RVassemblydir='K'.$RVk;
	$RVindex='MyOutput' unless (defined $RVindex and $RVindex=~/^\S+$/);
	if ($RVk>0 and defined $RVreadset) {
		$velveth_cmd=$path_velveth.' '.$RVassemblydir.' '.$RVk.' '.$RVreadset.' >> velveth.log';
		unless (&exec_cmd($velveth_cmd)) {
			print STDERR "Error: velveth failed: $RVindex $RVassemblydir\n";
			return 0;
		}
		$velvetg_cmd=$path_velvetg.' '.$RVassemblydir.' -exp_cov auto -cov_cutoff auto -amos_file yes'.' '.$RVvelvetg_set.' >> velvetg.log';
		unless (&exec_cmd($velvetg_cmd)) {
			print STDERR "Error: velvetg failed: $RVindex $RVassemblydir\n";
			return 0;
		}
	}
	else {
		print "Info: the assembly for read set ($RVreadset) at Kmer ($RVk) failed\n";
	}
	unless (chdir "./$RVassemblydir") {
		print STDERR "Velvet Error: can not cd dir ./$RVassemblydir：$!\n";
		return 0;
	}
	my $RVcurdir=getcwd;
#	my $RVassembly_fa=glob "$RVcurdir/contigs.fa";
	my $RVassembly_fa="$RVcurdir/contigs.fa";
	unless (-s $RVassembly_fa) {
		print STDERR "Error: Can not find Velvet output: $RVassembly_fa\n";
		return 0;
	}
	if ($mummer_each_K) {
		unless (&MumStats($input_fasta, $RVassembly_fa, "$RVindex.$RVassemblydir.AllContigs")) {
			print STDERR "Error: MUMmerplot running failed: AllContigs $RVindex $RVassemblydir\n";
			return 0;
		}
	}
	&ContigsStats($RVassembly_fa, 'longest_contig.fa', 0);
	unless (-s "$RVcurdir/longest_contig.fa") {
		print STDERR "Error: Can not find longest contigs file: $RVcurdir/longest_contig.fa";
		return 0;
	}
	if ($mummer_each_K) {
		unless (&MumStats($input_fasta, "$RVcurdir/longest_contig.fa", "$RVindex.$RVassemblydir.LongestContig")) {
			print STDERR "Error: MUMmerplot running failed: LongestContig $RVindex $RVassemblydir\n";
			return 0;
		}
	}
	unless (chdir "../") {
		print STDERR "Velvet Error: can not cd dir ../：$!\n";
		return 0;
	}
	
	return 1;
}



###MumStats plot fasta to reference
#&MumStats($ref, $fasta, $name)
sub MumStats {
	my ($MSref, $MSfa, $MSoutput)=@_;
	print "\n\nMummerplot map $MSfa to $MSref, output prefix: $MSoutput\n";
	my $mum_cmd="$path_nucmer -maxmatch -p=$MSoutput $MSref $MSfa";
	unless (&exec_cmd($mum_cmd)) {
		print STDERR "Error: failed to run nucmer\n";
		return 0;
	}
	$mum_cmd="$path_deltafilter -q $MSoutput.delta > $MSoutput.filter.q";
	unless (&exec_cmd($mum_cmd)) {
		print STDERR "Error: failed to run delta-filter\n";
		return 0;
	}
#	$mum_cmd="$path_mummerplot --large --layout -p $MSoutput --png $MSoutput.filter.q > /dev/null 2>&1";
	$mum_cmd="$path_mummerplot --large --layout -p $MSoutput --postscript $MSoutput.filter.q";
	unless (&exec_cmd($mum_cmd)) {
		print STDERR "Error: failed to run mummerplot\n";
		return 0;
	}
	
	return 1;
}



###Extract longest contig from contigs.fa file
#&ContigsStats(input, longest_output, Options)
#options: 0, 'sum', 'stats', 'n50' 
sub ContigsStats {
	my ($CSinput, $CSoutput, $CSoptions)=@_;
	my $CSseqid='';
	my $CScur_seq='';
	my $CSmax_id='';
	my $CSmax_seq='';
	my @CSseq_length=();
	open (CSIN, $CSinput) || die "Can not find input when extracting longest contig\n";
	while (my $CSline=<CSIN>) {
		chomp $CSline;
		if ($CSline=~/^>(\S+)\s*/) {
			if ($CSseqid ne '' and length($CSmax_seq)<length($CScur_seq)) {
				$CSmax_id=$CSseqid;
				$CSmax_seq=$CScur_seq;
			}
			push (@CSseq_length, length($CScur_seq));
			$CSseqid=$1;
			$CScur_seq='';
		}
		else {
			$CScur_seq.=$CSline;
		}
	}
	if ($CSseqid ne '' and length($CSmax_seq)<length($CScur_seq)) {
		$CSmax_id=$CSseqid;
		$CSmax_seq=$CScur_seq;
	}
	push (@CSseq_length, length($CScur_seq));
	@CSseq_length=sort {$b <=> $a} @CSseq_length;
	close CSIN;
	if ($CSoutput) {
		unlink($CSoutput) if (-e $CSoutput);
		open (CSOUT, ">>$CSoutput") || die "Can not write to file: $CSoutput\n";
		print CSOUT ">$CSmax_id\n";
		print CSOUT $CSmax_seq."\n";
		close CSOUT;
	}
	elsif ($CSoptions eq 'sum' or $CSoptions eq 'stats'  or $CSoptions eq 'n50') {
		my $CSsum=0;
		foreach (@CSseq_length) {
			$CSsum+=$_;
		}
		my $CSnum_seqs=scalar(@CSseq_length);
		if ($CSoptions eq 'sum') {
			return ($CSsum, $CSnum_seqs);###return total length and number of contigs
		}
		elsif ($CSoptions eq 'n50' or $CSoptions eq 'stats') {
			my $CSn50=&QuartileCalc(50, $CSsum, \@CSseq_length);
			if ($CSoptions eq 'n50') {
				return $CSn50;###return N50
			}
			elsif ($CSoptions eq 'stats') {
				my $CSnum_200=&CountNumSeqs(200, \@CSseq_length);###Number of sequences with length>=200
				my $CSnum_500=&CountNumSeqs(500, \@CSseq_length);###Number of sequences with length>=500
				my $CSnum_n50=&CountNumSeqs($CSn50, \@CSseq_length);###Number of sequences with length>=N50
				my $CSmin=min @CSseq_length;###Minimum sequence length
				my $CSn10=&QuartileCalc(10, $CSsum, \@CSseq_length);###N10
				my $CSn20=&QuartileCalc(20, $CSsum, \@CSseq_length);###N20
				my $CSn30=&QuartileCalc(30, $CSsum, \@CSseq_length);###N30
				my $CSn40=&QuartileCalc(40, $CSsum, \@CSseq_length);###N40
				my $CSn60=&QuartileCalc(60, $CSsum, \@CSseq_length);###N60
				my $CSn70=&QuartileCalc(70, $CSsum, \@CSseq_length);###N70
				my $CSn80=&QuartileCalc(80, $CSsum, \@CSseq_length);###N80
				my $CSn90=&QuartileCalc(90, $CSsum, \@CSseq_length);###N90
				my $CSmax=length($CSmax_seq);###Maximum sequence length
				return ($CSsum, $CSnum_seqs, $CSnum_200, $CSnum_500, $CSnum_n50, $CSn10, $CSn20, $CSn30, $CSn40, $CSn50, $CSn60, $CSn70, $CSn80, $CSn90, $CSmax, $CSmax_id, $CSinput);
			}
		}
	}
}
###calculate N50
sub QuartileCalc {
	my ($QCquartile, $QCsum, $QCnums) = @_;
	unless ($QCquartile>=0 and $QCquartile<=100 and $QCquartile<=$QCsum and scalar(@{$QCnums})>1) {
		die "Error in calculating contigs Quartiles\n";
	}
	my $QCcount=0;
	foreach my $QCseq_length_ind (@{$QCnums}) {
		$QCcount+=$QCseq_length_ind;
		if ($QCcount >= $QCsum*$QCquartile/100){
			return ($QCseq_length_ind);
			last;
		}
	}
}
###CountSeqs
#&CountSeqs($min, @array)
sub CountNumSeqs {
	my ($CNSquartile, $CNSlengths)= @_;
	my $CNScount=0;
	foreach my $CNSseq_length_ind (@{$CNSlengths}) {
		if ($CNSseq_length_ind >=$CNSquartile){
			$CNScount++;
		}
	}
	return ($CNScount);
}



###Step3 summarize the contigs
#&AssemblySummary()
#global variables: 
sub AssemblySummary {
	
	local *ASSTATS;

	unless (chdir "$output_dir") {
		print STDERR "Error: Summary error: can not cd dir $output_dir：$!\n";
	}
	my $newdir=$output_dir.'/3.summary';
	if (-d $newdir) {
		unless (&exec_cmd("rm -rf $newdir > /dev/null 2>&1")) {
			print STDERR "Warnings: failed to delete dir: $newdir\n";
		}
	}
	unless (mkdir ($newdir, 0766)) {
		print STDERR "Error: Summary error: can not create directory $newdir\n";
		return 0;
	}
	unless (chdir "$newdir") {
		print STDERR "Error: Summary error: can not cd dir $newdir：$!\n";
		return 0;
	}
	my @AScontig_files=glob "$output_dir/2.assembly/*/*contigs.fa";
	unless (scalar(@AScontig_files)>0) {
		print STDERR "Error: Summary error: Can not find $output_dir/2.assembly/*/*contigs.fa";
		return 0;
	}
	unlink ("FinalStats.txt") if (-e "FinalStats.txt");
	unless (open (ASSTATS, ">>FinalStats.txt")) {
		print STDERR "Error: Summary error: can not write to file: $output_dir'/3.summary/FinalStats.txt'";
		return 0;
	}
	foreach my $ASindex (keys %RCFcomb) {
		my $ASnum_libs=scalar(@{$RCFcomb{$ASindex}});
		print ASSTATS "Index $ASindex Accepted: \n";
		for (my $i=0; $i<$ASnum_libs; $i++) {
			print ASSTATS "\t@{${$RCFcomb{$ASindex}}[$i]}\n";
		}
	}
	print "\nNum_bp\tNum_seqs\tNum>200\tNum>500\tNum>n50\tn10\tn20\tn30\tn40\tn50\tn60\tn70\tn80\tn90\tmax\tmax_id\tFile\n";
	foreach my $AScontigs (@AScontig_files) {
		print ASSTATS join("\t", &ContigsStats($AScontigs, 0, 'stats')),"\n";
	}
	close ASSTATS;
	
	unless (&exec_cmd("cp -f $output_dir/2.assembly/*/*.png $output_dir/2.assembly/*/*.ps $newdir/")) {
		print STDERR "Warnings: failed to collect images\n";
	}
	if (mkdir ("$newdir/all", 0766)) {
		foreach (@allcontigfile) {
			unless (&exec_cmd("cp -f $_ $newdir/all/")) {
				print STDERR "Warnings: copy file failed: $_ to $newdir/all/\n";
			}
		}
	}
	else {
		print STDERR "Warnings: Summary error: can not create directory $newdir\n";
	}
	
	if (mkdir ("$newdir/longest", 0766)) {
		foreach (@longestcontigfiles) {
			unless (&exec_cmd("cp -f $_ $newdir/longest/")) {
				print STDERR "Warnings: copy file failed: $_ to $newdir/longest/\n";
			}
		}
	}
	else {
		print STDERR "Warnings: Summary error: can not create directory $newdir\n";
	}
	
	unless (chdir "$output_dir") {
		print STDERR "Warnings: Summary error: can not cd dir $output_dir：$!\n";
	}
}



sub FastaSum {
	my $FSfile=shift @_;
	
	local *FSINPUT;
	unless (open (FSINPUT, "< $FSfile")) {
		print STDERR "Error: Can not open fasta files specified: $FSfile\n";
		return 0;
	}
	my $FSnum_seq=0;
	my $FStotal_length=0;
	my $FSseqid='';
	my $FSline_num=0;
	while (my $FSline=<FSINPUT>) {
		chomp $FSline;
		my $FSline_num++;
		if ($FSline =~/^>/) {
			$FSnum_seq++;
			$FSseqid=$FSline;
		}
		else {
			$FStotal_length+=length($FSline);
			if ($FSline=~/[^atcgATCGnN]/) {
				print STDERR "Warnings: Non-ATCGatcgNn letters at line $FSline_num in sequence $FSseqid in file $FSfile\n";
			}
		}
	}
	close FSINPUT;
#	return ($FSnum_seq, $FStotal_length);
	return ($FStotal_length);
}



###Get (path and extesion name)-removed fileID
#&RetrvNoExt(file)
sub RetrvNoExt {
	my $RNE_ori=shift @_;
	chomp $RNE_ori;
	my ($RNE_new, $RNE_base)=('', '');
	($RNE_base=$RNE_ori)=~ s/.*\///s; 
	($RNE_new=$RNE_base)=~s/^(.*)\.\w+$/$1/;
	return $RNE_new;
}



sub mytime {
	my($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime();
	$year += 1900;
	$mon  += 1;
	my $btime = sprintf("%04d%02d%02d %02d:%02d:%02d",$year,$mon,$day,$hour,$min,$sec);
	return $btime;
}



sub exec_cmd {
	my ($cmd) = @_;
	print "\#" x 70 ."\n".&mytime()."CMD: $cmd\n";
	my $start_time = time();
	my $return_code = system($cmd);
	my $end_time = time();
#	if ($return_code == -1) {
#		print “failed to execute: $!\n”;
#	}
#	elsif ($return_code & 127) {
#		printf “child died with signal %d, %s coredump\n”, ($return_code & 127),  ($return_code & 128) ? ‘with’ : ‘without’;
#	}
#	else {
#		printf “child exited with value %d\n”, $return_code >> 8;
#	}
	if ($return_code) {
		print "Error, cmd: $cmd died with ReturnCode $return_code\n";
		return 0;
	}
	else {
		print "\nFinished command: $cmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n\n";
		return 1;
	}
}



###delete path
#my $curDir=getcwd;
#&DeletePath(PATH);
#chdir $curDir || die "Error when deleting directory: $curDir\n";
sub DeletePath{
	my $DPpath = shift @_;
	chdir $DPpath || die "Error when deleting directory: $DPpath\n";
	#get all the files in that directory.
	my @DPfiles=<*>;
	foreach(@DPfiles){
		if(-d $_){
		#if the destination file is a directory, go recursion.
		&DeletePath($_);
	}else{
		unlink;
	}
	}
	#Go up and del the destination directory.
	chdir "../";
	rmdir $DPpath;
}
