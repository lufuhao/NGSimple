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
	LUFUHAO20140915

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



our ($help_message, $verbose, $version);
our ($config_file, $output_dir, $test_run_stepwiseK, $num_threads);
our ($path_java, $path_trimmomatic, $trimmo_minlen);
our ($input_fasta, $path_velvetg, $path_velveth);
our ($path_nucmer, $path_deltafilter, $path_mummerplot);

GetOptions (
'help|h!'=> \$help_message,
'config|c:s' => \$config_file,
'reference|r:s' => \$input_fasta,
'threads|t:i' => \$num_threads,
'output_dir|d:s' => \$output_dir,
'path_java:s' => \$path_java,
'path_trimmomatic:s' => \$path_trimmomatic,
'minlength|l:i' => \$trimmo_minlen,
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
our (@lib_names, @lib_types, @lib_read_lengths, @lib_covs, @lib_sizes, %RCFlib, %RCFcomb, $MaxVelvetCategories);
$num_threads=1 unless (defined $num_threads);
our $experence_factor=0.2;
our $cmd='';
$path_java='java' unless (defined $path_java);
$path_trimmomatic="$Bin/utils/Trimmomatic/v0.32/x86_64/trimmomatic-0.32.jar" unless (defined $path_trimmomatic);
our $quality_score_format= 'phred33' unless (defined $quality_score_format);
die "Quality format accept \'phred33\' or \'phred64\' only\n" unless ($quality_score_format eq 'phred33' or $quality_score_format eq 'phred64');
$test_run_stepwiseK=1 unless (defined $test_run_stepwiseK);
$trimmo_minlen=36;
$MaxVelvetCategories=10;
$path_velveth='velveth' unless (defined $path_velveth);
$path_velvetg='velvetg' unless (defined $path_velvetg);
$path_nucmer='nucmer' unless (defined $path_nucmer);
$path_deltafilter='delta-filter' unless (defined $path_deltafilter);
$path_mummerplot='mummerplot' unless (defined $path_mummerplot);
###test output
our $run_dir=getcwd;
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
foreach my $AssemblyIndex (keys %RCFcomb) {
	print "Assembling Index: $AssemblyIndex\n";
	&SettingVelvet($AssemblyIndex);
}
&AssemblySummary();
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
}



###Step2: Set velvet
#&SetVelvet($index)
#global variables: %RCFcomb, $test_run_stepwiseK, $input_fasta, $experence_factor
sub SettingVelvet {
	my $SVindex=shift @_;
	chdir "$output_dir" or die "Velvet error: can not cd dir $output_dir：$!\n";
	my $newdir="$output_dir/2.assembly/$SVindex";
	&exec_cmd("rm -rf $newdir") if (-d $newdir);
	mkdir ($newdir, 0766) || die "Velvet error: can not create directory $newdir\n";
	chdir "$newdir" or die "Vetvet error: can not cd dir $newdir：$!\n";
	my @SVseq_files=();
	my $SVreadset='';
	my $SV_num_lib=1;
	my $SVtotal_coverage=0;
	my $SVmin_readlength=10000;
	for (my $SVi=0; $SVi<scalar(@{$RCFcomb{$SVindex}}); $SVi++) {
		if (${${$RCFcomb{$SVindex}}[$SVi]}[0] eq 'se') {
			my $SVfastq_file='';
			($SVfastq_file)=&ReadFastqFiles(${${$RCFcomb{$SVindex}}[$SVi]}[0], ${${$RCFcomb{$SVindex}}[$SVi]}[3], ${${$RCFcomb{$SVindex}}[$SVi]}[2], ${${$RCFcomb{$SVindex}}[$SVi]}[1]);
			push (@SVseq_files, $SVfastq_file);
			$SVreadset.=" -short -fastq $SVfastq_file";
		}
		elsif (${${$RCFcomb{$SVindex}}[$SVi]}[0] eq 'pe' or ${${$RCFcomb{$SVindex}}[$SVi]}[0] eq 'mp') {
			my ($SVfastq_R1, $SVfastq_R2)=('', '');
			($SVfastq_R1, $SVfastq_R2)=&ReadFastqFiles(${${$RCFcomb{$SVindex}}[$SVi]}[0], ${${$RCFcomb{$SVindex}}[$SVi]}[3], ${${$RCFcomb{$SVindex}}[$SVi]}[2], ${${$RCFcomb{$SVindex}}[$SVi]}[1]);
			push (@SVseq_files, $SVfastq_R1);
			push (@SVseq_files, $SVfastq_R2);
			if ($SV_num_lib==1) {
				$SVreadset.=" -shortPaired -fastq -separate $SVfastq_R1 $SVfastq_R2 ";
			}
			else {
				$SVreadset.=" -shortPaired$SV_num_lib -fastq -separate $SVfastq_R1 $SVfastq_R2 ";
			}
			$SV_num_lib++;
		}
		else {
			die "Unknown libary type: ${${$RCFcomb{$SVindex}}[$SVi]}[0]\n";
		}
		$SVmin_readlength=min $SVmin_readlength,${${$RCFcomb{$SVindex}}[$SVi]}[1];
		$SVtotal_coverage+=${${$RCFcomb{$SVindex}}[$SVi]}[2];
	}
	my $SVgenome_size=&fastasum($input_fasta);
	my $SVtargetKmer=$SVtotal_coverage*$experence_factor;
	print "\n\nReference size $SVgenome_size\nExpect K-mer:$SVtargetKmer\nAssembly fies: @SVseq_files\n";
	my $SVguess_beskK=&GuessBestKmer($SVgenome_size, $SVtargetKmer, @SVseq_files);
	print "Guess bast K-mer: $SVguess_beskK\n\n";
	chdir "$newdir" or die "Vetvet error: can not cd dir $newdir：$!\n";
	if ($test_run_stepwiseK==0) {
		&RunVelvet($SVguess_beskK, $SVreadset);
		link("./K$SVguess_beskK/contigs.fa", "$SVindex.K$SVguess_beskK.contigs.fa");
		link("./K$SVguess_beskK/longest_contig.fa", "$SVindex.K$SVguess_beskK.LongestContig.fa");
		link("./K$SVguess_beskK/AllContigs.png", "$SVindex.K$SVguess_beskK.AllContigs.png");
		link("./K$SVguess_beskK/LongestContig.png", "$SVindex.K$SVguess_beskK.LongestContig.png");
	}
	else {
		my $SVrange=int($SVguess_beskK/5);
		$SVrange+=1 if ($SVrange%2==1);
		my $SVbestKmer=&SelectBestAssembly($SVguess_beskK-$SVrange, $SVguess_beskK+$SVrange, 4, $SVreadset);
		my $SVbestKmer2=&SelectBestAssembly($SVbestKmer-2, $SVbestKmer=+2, 2, $SVreadset);
		link("./K$SVbestKmer2/contigs.fa", "$SVindex.K$SVbestKmer2.contigs.fa");
		link("./K$SVbestKmer2/longest_contig.fa", "$SVindex.K$SVbestKmer2.LongestContig.fa");
		link("./K$SVbestKmer2/AllContigs.png", "$SVindex.K$SVbestKmer2.AllContigs.png");
		link("./K$SVbestKmer2/LongestContig.png", "$SVindex.K$SVbestKmer2.LongestContig.png");
	
	}
	chdir "$output_dir" or die "Velvet error: can not cd dir $output_dir：$!\n";
}



###Step3 summarize the contigs
#&AssemblySummary()
#global variables: 
sub AssemblySummary {
	chdir "$output_dir" or die "Summary error: can not cd dir $output_dir：$!\n";
	my $newdir=$output_dir.'/3.summary';
	&exec_cmd("rm -rf $newdir") if (-d $newdir);
	mkdir ($newdir, 0766) || die "Summary error: can not create directory $newdir\n";
	chdir "$newdir" or die "Summary error: can not cd dir $newdir：$!\n";
	my @AScontig_files=glob "$output_dir/2.assembly/*/*contigs.fa";
	die "Summary error: Can not find $output_dir/2.assembly/*/*contigs.fa" unless (@AScontig_files);
	open (STATS, ">>FinalStats.txt") || die "Summary error: can not write to file: $output_dir'/3.summary/FinalStats.txt'";
	foreach my $ASindex (keys %RCFcomb) {
		my $ASnum_libs=scalar(@{$RCFcomb{$ASindex}});
		print STATS "Index $ASindex Accepted: \n";
		for (my $i=0; $i<$ASnum_libs; $i++) {
			print STATS "\t@{${$RCFcomb{$ASindex}}[$i]}\n";
		}
	}
	
	foreach my $AScontigs (@AScontig_files) {
		print STATS &ContigsStats($AScontigs, 0, 'stats'), "\n";
	}
	&exec_cmd("cp $output_dir/2.assembly/*/*.png ./");
	chdir "$output_dir" or die "Summary error: can not cd dir $output_dir：$!\n";
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
				die "Can not find Fq files: $RFFfile_basename'_1.fq' and $RFFfile_basename'_2.fq'";
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
	if (lc($quality_score_format) eq 'phred33') {
		$cmd.=' -phred33 ';
	}
	elsif (lc($quality_score_format) eq 'phred64') {
		$cmd.=' -phred64 ';
	}
	else {
		die "Trimmomatic error: Uknown quality score format $quality_score_format\n";
	}
	$cmd.=$TMfile_setting."HEADCROP:8 LEADING:18 TRAILING:18 SLIDINGWINDOW:4:18 MINLEN:$trimmo_minlen";
	&exec_cmd($cmd);
	chdir "$output_dir" or die "Trimmomatic error: can not cd dir $output_dir：$!\n";
	print "###Trimmomatic END###\n\n\n";
	if ($TMfq2_output ne '' and -s $TMfq2_output) {
		return ($TMfq1_output, $TMfq2_output);
	} elsif ($TMfq2_output eq '') {
		return ($TMfq1_output);
	}
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



###RunVelvet
#&RunVelvet($k, $readset)
#global variables: $path_velveth, $path_velvetg
sub RunVelvet {
	my ($RVk, $RVreadset)=@_;
	my $velveth_cmd='';
	my $velvetg_cmd='';
	my $RVassemblydir='K'.$RVk;
	if ($RVk>0 and defined $RVreadset) {
		$velveth_cmd=$path_velveth.' '.$RVassemblydir.' '.$RVk.' '.$RVreadset;
		&exec_cmd($velveth_cmd);
		$velvetg_cmd=$path_velvetg.' '.$RVassemblydir.' -exp_cov auto -cov_cutoff auto -amos_file yes';
		&exec_cmd($velvetg_cmd);
	}
	else {
		print "The assembly for read set ($RVreadset) at Kmer ($RVk) failed\n";
	}
	chdir "./$RVassemblydir" || die "Velvet Error: can not cd dir ./$RVassemblydir：$!\n";
	my $RVassembly_fa=glob "contigs.fa";
	&MumStats($input_fasta, $RVassembly_fa, 'AllContigs');
	&ContigsStats($input_fasta, 'longest_contig.fa', 0);
	&MumStats($input_fasta, 'longest_contig.fa', 'LongestContig');
	chdir "../" || die "Velvet Error: can not cd dir ../：$!\n";
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
		print CSOUT '>'.$CSmax_id."\n".$CSmax_seq."\n";
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
			my $CSn50=&QuartileCalc(50, $CSsum, @CSseq_length);
			if ($CSoptions eq 'n50') {
				return $CSn50;###return N50
			}
			elsif ($CSoptions eq 'stats') {
				my $CSnum_200=&CountNumSeqs(200, @CSseq_length);###Number of sequences with length>=200
				my $CSnum_500=&CountNumSeqs(500, @CSseq_length);###Number of sequences with length>=500
				my $CSnum_n50=&CountNumSeqs($CSn50, @CSseq_length);###Number of sequences with length>=N50
				my $CSmin=min @CSseq_length;###Minimum sequence length
				my $CSn10=&QuartileCalc(10, $CSsum, @CSseq_length);###N10
				my $CSn20=&QuartileCalc(20, $CSsum, @CSseq_length);###N20
				my $CSn30=&QuartileCalc(30, $CSsum, @CSseq_length);###N30
				my $CSn40=&QuartileCalc(40, $CSsum, @CSseq_length);###N40
				my $CSn60=&QuartileCalc(60, $CSsum, @CSseq_length);###N60
				my $CSn70=&QuartileCalc(70, $CSsum, @CSseq_length);###N70
				my $CSn80=&QuartileCalc(80, $CSsum, @CSseq_length);###N80
				my $CSn90=&QuartileCalc(90, $CSsum, @CSseq_length);###N90
				my $CSmax=length($CSmax_seq);###Maximum sequence length
				return ($CSsum, $CSnum_seqs, $CSnum_200, $CSnum_500, $CSnum_n50, $CSn10, $CSn20, $CSn30, $CSn40, $CSn50, $CSn60, $CSn70, $CSn80, $CSn90, $CSmax, $CSmax_id, $CSinput);
			}
		}
	}
}
###calculate N50
sub QuartileCalc {
	my ($QCquartile, $QCsum, @QCnums) = @_;
	die "Error in calculating contigs Quartiles\n " unless ($QCquartile>=0 and $QCquartile<=100 and $QCquartile<=$QCsum and scalar(@QCnums)>1);
	my $QCcount=0;
	foreach my $QCseq_length_ind (@QCnums) {
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
	my ($CNSquartile, @CNSlengths)= @_;
	my $CNScount=0;
	foreach my $CNSseq_length_ind (@CNSlengths) {
		if ($CNSseq_length_ind >=$CNSquartile){
			$CNScount++;
		}
	}
	return ($CNScount);
}


###select largest n50 for ./K*/contigs.fa
#&SelectBestAssembly(min, max, step, $SVreadset)
sub SelectBestAssembly {
	my ($SBAminK, $SBAmaxK, $SBAstep, $SBAreadset)=@_;
	for (my $SBAj=$SBAminK; $SBAj<=$SBAmaxK; $SBAj+=$SBAstep) {
		&RunVelvet($SBAj, $SBAreadset) unless (-d "./K$SBAj");
	}
	my $SBAbestK=0;
	my $SBAmax_n50=0;
	my @SBAfiles=glob "./K*/contigs.fa";
	die "Can not find any ./K*/contigs.fa velvet assemblies" unless (scalar(@SBAfiles)>0);
	foreach my $SBAcontig_file (@SBAfiles) {
		my $SBAind_n50=&ContigsStats($SBAcontig_file, 0, 'n50');
		if ($SBAind_n50>$SBAmax_n50) {
			if ($SBAcontig_file=~m/\/\S+\/K(\d{1,3})\//) {
				$SBAbestK=$1;
			}else { die "Can not find best K for velvet\n";}
		}
	}
	if (defined $SBAbestK and $SBAbestK >0) {
		return ($SBAbestK);
	}else {
		die "Can not calcaulate best Kmer between $SBAminK - $SBAmaxK with a step $SBAstep dor read set \'$SBAreadset\'\n";
	}
	
}



###MumStats plot fasta to reference
#&MumStats($ref, $fasta, $name)
sub MumStats {
	my ($MSref, $MSfa, $MSoutput)=@_;
	my $mum_cmd="$path_nucmer -maxmatch -p=$MSoutput $MSref $MSfa";
	&exec_cmd($mum_cmd);
	$mum_cmd="$path_deltafilter -q $MSoutput.delta > $MSoutput.filter.q";
	&exec_cmd($mum_cmd);
	$mum_cmd="$path_mummerplot --large --layout -p=$MSoutput --png $MSoutput.filter.q ";
	&exec_cmd($mum_cmd);
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
	die "Can not find any read\n" if (scalar(@GBKread_length)<1);
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
#	return ($FSnum_seq, $FStotal_length);
	return ($FStotal_length);
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
	print &mytime()."CMD: $cmd\n";
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
		return $return_code;
	}
	else {
		print "\nFinished command: $cmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n\n";
		return $return_code;
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
