#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]

Version: 
	v20141119

Requirements:
	Programs:
	Modiles: Scalar::Util, Cwd, Getopt::Long, FindBin

Descriptions:
	Determine the insert size given pairs of seqing data by
	mapping them to a reference.

Options:
	--help|-h
		Print this help/usage;
	--FastqR1|-1	<Fastq/Fastq_list>
		[Optional] Input read1 list separated by comma;
	--FastqR2|-2	<File/Fastq_list>
		[Optional] Input read2 list separated by comma;
	--bam|-b <Bam/Bam_List>
		[Optional] Input bam list separated by comma;
	--procedure|-p	<list>
		[Necessary] procedure list separated by comma;
		1.cutadapt #Remove adatptors from 5' or 3' end;
		2.trimmomatic #trimming and filtering
		3.fastq-join #Detect if 2 reads in one pair overlaps;
		4.fastx_trimmer #Trim (chimeric) reads to certain length;
		5.fastx_reverse_complement #rc RF to FR;
		6.bwa #map reads to reference;
		7.Collect and report insert size;
		Defualt: 1,2,6,7;
	--fqformat <phred33|phred64>
		Fastq quality score format; Default: phred33;
	--threads|-t	<INT>
	--prefix	<STR>
		Prefix of output file name; Default: 'Group';
	--path_java	</path/to/java>
    #fastqc
	--run_fastqc
		Whether to run fastqc for 1,2,3,4,5; Default: No;
	--path_fastqc	</path/to/fastqc>
    #cutadapt
	--path_cutadapt	</path/to/cutadapt>
	--cutadapt_overlap	<INT>
		Minimum overlap length for --overlap;
	--cutadapt_times	<INT>
		Try to remove adapters at most INT times for --times;
	--file_adaptor	</path/to/adators_file>
    #trimmomatic
	--path_trimmomatic:s	</path/to/trimmomatic.jar>
	--contaminants	</path/to/seqprimer_file>
	--trimmomatic_minlen	<INT>
		Minimum read legnth to keep a read when trimming;
    #fastqjoin
	--path_fastqjoin	</path/to/fastq-join>
	--min_overlap	<INT>
		Minimum overlap length when joining two read;
	--sizebin	<INT>
		Size bin for distribute the insert sizes;
    #trimmer
	--path_fastxtrimmer	</path/to/fastx_trimmer>
	--trimmer_minlen	<INT>
		Trim all reads to INT bp length;
    #reversecomplement
	--path_fastxrc	</path/to/fastx_reverse_complement>
    #bwa
	--path_bwa	</path/to/bwa>
	--reference|-r	</path/to/reference_fasta>
	--minimum_insert	<INT>
		Minimum_insert size for paired reads;
		Defaults: 0;
	--maximum_insert	<INT>
		Maximum_insert size for paired reads;
	--min_mapq	<INT>
		Minimum mapping quality value to keep;
	--reindex
		Delete old index, and Re-run bwa index;
    #samtools
	--path_samtools	</path/to/samtools>
    #picard
	--path_picard	</path/to/PICARD_HOME>
	--picard_memory	<STR>
		Maximum memory for java -Xmx<STR>, Default: 2g;
	--verbose
		Detailed output for trouble-shooting;
	--version|v!
		Print current SCRIPT version;

Example:
	perl $0 -1 forward.fq -2 reverse.fq --procedure 2,6,7 \
	--fqformat phred33 \
	--contaminants primers.fa \
	--reference reference.fa \
	--maximum_insert 60000
	

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk

EOH
###HELP ends##########################
die USAGE unless @ARGV;



###Receving parameter#################
our ($help, $FastqR1_list, $FastqR2_list, $bamlist, $procedure, $quality_score_format, $num_threads, $prefix, $path_java, $verbose, $ver);
our ($test_run_fastqc, $path_fastqc);
our ($path_cutadapt, $cutadapt_overlap, $cutadapt_times, $file_adaptor);
our ($path_trimmomatic, $contaminants, $trimmo_minlen);
our ($path_fastqjoin, $min_overlap, $sizebin);
our ($path_fastxtrimmer, $trimmer_minlen);
our ($path_fastxreversecomplement);
our ($path_bwa, $bwaref, $minimum_insertsize, $maximum_insertsize, $min_mapq, $test_reIndex);
our ($path_samtools);
our ($path_picard, $picard_memory);

GetOptions(
	"help|h!" => \$help,
	"FastqR1|1:s" => \$FastqR1_list,
	"FastqR2|2:s" => \$FastqR2_list,
	"bam|b:s" => \$bamlist,
	"procedure|p:s" => \$procedure,
	"fqformat:s" => \$quality_score_format,
	"threads|t:i" => \$num_threads,
	"prefix:s" => \$prefix,
	"path_java:s" => \$path_java,
###fastqc
	"path_fastqc:s" => \$path_fastqc,
	"run_fastqc!" => \$test_run_fastqc,
###cutadapt
	"path_cutadapt:s" => \$path_cutadapt,
	"cutadapt_overlap:i" => \$cutadapt_overlap,
	"cutadapt_times:i" => \$cutadapt_times,
	"file_adaptor:s" => \$file_adaptor,
###trimmomatic
	"path_trimmomatic:s" => \$path_trimmomatic,
	"contaminants:s" => \$contaminants,
	"trimmomatic_minlen:i" => \$trimmo_minlen,
###fastqjoin
	"path_fastqjoin:s" => \$path_fastqjoin,
	"min_overlap:i" => \$min_overlap,
	"sizebin:i" => \$sizebin,
###trimmer
	"path_fastxtrimmer:s" => \$path_fastxtrimmer,
	"trimmer_minlen:i" => \$trimmer_minlen,
###reversecomplement
	"path_fastxrc:s" => \$path_fastxreversecomplement,
###bwa
	"path_bwa:s" => \$path_bwa,
	"reference|r:s" => \$bwaref,
	"minimum_insert:i" => \$minimum_insertsize,
	"maximum_insert:i" => \$maximum_insertsize,
	"min_mapq:i" => \$min_mapq,
	"reindex!" => \$test_reIndex,
###$path_samtools
	"path_samtools:s" => \$path_samtools,
###picard
	"path_picard:s" => \$path_picard,
	"picard_memory:s" => \$picard_memory,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;



###Default##############################################
our $cmd='';
our $stepth=0;
our $Stage='Stage '.$stepth.' Loading';

print "\n\n\n###BEGIN###$Stage\n";
$prefix='Group' unless (defined $prefix);
our $input_format='fastq' unless (defined $input_format);

#defined rundir
our $pwd_dir=getcwd;
our $run_dir=$pwd_dir unless (defined $run_dir);
unless (-d $run_dir) {
	mkdir ($run_dir, 0766) || die "Error at $Stage: can not create path: $run_dir\n";
}
print "Current working path was set to $run_dir\n" if (defined $verbose);

#Fastq_input_check
our $test_use_reads=0;
our $test_use_bam=0;
our @ReadGroups=();
our $num_groupofreads=0;
if (defined $FastqR1_list and defined $FastqR2_list) {
	our @FastqR1_arr=split(/,/, $FastqR1_list);
	our @FastqR2_arr=split(/,/, $FastqR2_list);
	@FastqR1_arr=&AddFilePath(@FastqR1_arr);
	@FastqR2_arr=&AddFilePath(@FastqR2_arr);
	if (scalar(@FastqR1_arr) == scalar(@FastqR2_arr)) {
		$num_groupofreads=scalar(@FastqR1_arr);
		unless ($num_groupofreads>0) {
			die "Error at $Stage: Please specify input of paired reads\n";
		}
	}
	else {
		die "Error at $Stage: Un-equal Paired number for R1 and R2\n";
	}
	for (my $i=0; $i<$num_groupofreads;$i++) {
		unless (-s $FastqR1_arr[$i]) {
			die "Error at $Stage: The ",$i+1,"th file in R1 does not exist\n";
		}
		unless (-s $FastqR2_arr[$i]) {
			die "Error at $Stage: The ",$i+1,"th file in R2 does not exist\n";
		}
		&AddToReadGroup($i, $FastqR1_arr[$i], $FastqR2_arr[$i], 0);
		$test_use_reads=1;
	}
	map {print "Group: \n1. ".$_->[0]."\n2. ".$_->[1]."\n"} @ReadGroups;
}


#Bam_input_check
our @BamGroups=();
if ($test_use_reads==0) {
	if (defined $bamlist) {
		@BamGroups=split(/,/,$bamlist);
		@BamGroups=&AddFilePath(@BamGroups);
		$num_groupofreads=scalar(@BamGroups);
		unless ($num_groupofreads>0){
			die "Error at $Stage: Please specify input of paired reads\n";
		}
		foreach (@BamGroups) {
			unless (-s $_) {
				die "Error at $Stage: Some file in bamlist does not exist\n";
			}
		}
		$test_use_bam=1;
	}
}

###Parameter design for each procedures##############################
our $test_run_cutadapt=0;
our $test_run_trimmomatic=0;
our $test_run_fastqjoin=0;
our $test_run_trimmer=0;
our $test_run_reversecomplement=0;
our $test_run_bwa=0;
#our $test_run_extractpaired=0;
#our $test_run_extractunique=0;
our $test_run_sizecollect=0;
$procedure="1,2,6,7" unless (defined $procedure);
#print "Procedure: $procedure\n";                      ###test
our @procdures=split(/,/,$procedure);
foreach (@procdures) {
	MYTEST02: {
		if ($_==1) {$test_run_cutadapt++; last MYTEST02};
		if ($_==2) {$test_run_trimmomatic++; last MYTEST02};
		if ($_==3) {$test_run_fastqjoin++; last MYTEST02};
		if ($_==4) {$test_run_trimmer++; last MYTEST02};
		if ($_==5) {$test_run_reversecomplement++; last MYTEST02};
		if ($_==6) {$test_run_bwa++; last MYTEST02};
#		if ($_==7) {$test_run_extractpaired++; last MYTEST02};
#		if ($_==8) {$test_run_extractunique++; last MYTEST02};
		if ($_==7) {$test_run_sizecollect=1; last MYTEST02};
		die "Error at $Stage: Unknown Procedure parameter: $_\n";
	}
}
if ($test_use_bam) {
	$test_run_cutadapt=0;
	$test_run_trimmomatic=0;
	$test_run_fastqjoin=0;
	$test_run_trimmer=0;
	$test_run_reversecomplement=0;
	$test_run_bwa=0
}

print "Checking Programs and their parameters...\n";
###fastqc
$test_run_fastqc=0 unless (defined $test_run_fastqc);
if ($test_run_fastqc>0) {
	if (defined $path_fastqc){
		&TestCmdExist($path_fastqc);
	}
	else {
		$path_fastqc='fastqc';
	}
}

###cutadapt
if ($test_run_cutadapt) {
	defined $path_cutadapt ? (&TestCmdExist($path_cutadapt)) : ($path_cutadapt='cutadapt');
	unless (-s $file_adaptor) {
		die "Error at $Stage: Plase specify fasta file containing adaptors\n";
	}
	$cutadapt_overlap=8 unless (defined $cutadapt_overlap);
}

#test java
#if ($test_run_trimmomatic or $test_run_bwa or $test_run_extractpaired or $test_run_extractunique) {
if ($test_run_trimmomatic or $test_run_bwa){
	(defined $path_java) ? (&TestCmdExist($path_java)) : ($path_java='java');
}
#test trimmomatic path
if ($test_run_trimmomatic) {
#	die "Please specify Path/To/Trimmomatic file\n" unless (defined $path_trimmomatic);
	$path_trimmomatic="$Bin/utils/Trimmomatic/v0.32/x86_64/trimmomatic-0.32.jar" unless (defined $path_trimmomatic);
	unless ( -s $contaminants) {
		die "Error at $Stage: Can not find Trimmomatic contaminants files\n";
	}
	$quality_score_format='phred33' unless (defined $quality_score_format);
	$trimmo_minlen=36 unless (defined $trimmo_minlen);
}

#test PICARD path & samtools
#if ($test_run_bwa or $test_run_extractpaired or $test_run_extractunique) {
if ($test_run_bwa){
	unless (defined $path_picard) {
		if ($ENV{'PICARD_HOME'}) {
			$path_picard=$ENV{'PICARD_HOME'};
		}
		else {
			$path_picard="$Bin/utils/picard_tools/v1.108/picard-tools-1.108";
		}
	}
	unless (-e "$path_picard/MarkDuplicates.jar") {
		die "Error at $Stage: Can not find PICARD_HOME or path to PICARD\n";
	}
	$picard_memory='2g' unless (defined $picard_memory);
}
(defined $path_samtools) ? (&TestCmdExist($path_samtools)) : ($path_samtools='samtools');

###fastqjoin
if ($test_run_fastqjoin) {
	(defined $path_fastqjoin) ? (&TestCmdExist($path_fastqjoin)) : ($path_fastqjoin='fastq-join');
	$min_overlap=20 unless (defined $min_overlap);
	$sizebin=1000 unless (defined $sizebin);
}

#test fastx_reverse_complement
if ($test_run_reversecomplement) {
	(defined $path_fastxreversecomplement) ? (&TestCmdExist($path_fastxreversecomplement)) : ($path_fastxreversecomplement='fastx_reverse_complement');
}


#test fastx_trimmer
if ($test_run_trimmer) {
	(defined $path_fastxtrimmer) ? (&TestCmdExist($path_fastxtrimmer)) : ($path_fastxtrimmer='fastx_trimmer');
	die "Error at $Stage: Need to specified trimmer_minlen fir Fastx_Trimmer\n" unless (defined $trimmer_minlen);
}

###test bwa
if ($test_run_bwa) {
	(defined $path_bwa) ? (&TestCmdExist($path_bwa)) : ($path_bwa='bwa');
	unless (defined $bwaref and -s $bwaref) {
		die "Error at $Stage: can not find Fasta reference for BWA\n";
	}
	unless (defined $maximum_insertsize) {
		die "Error at $Stage: Please specify maximum insert size for BWA\n";
	}
	
	$test_reIndex=0 unless (defined $test_reIndex);
	die "Please specify maximum insert size for BWA\n" unless (defined $maximum_insertsize);
}
$min_mapq=2 unless (defined $min_mapq);
$minimum_insertsize=0 unless (defined $minimum_insertsize);



###Main#################################################
&PreCheck() unless ($test_use_bam);
foreach (@procdures){
	MYTEST01: {
		if ($_==1 and $test_use_bam==0){&CutAdapt(); last MYTEST01};
		if ($_==2 and $test_use_bam==0){&Trimmomatic(); last MYTEST01};
		if ($_==3 and $test_use_bam==0){&FastqJoin(); last MYTEST01};
		if ($_==4 and $test_use_bam==0){&FastxTrimmer(); last MYTEST01};
		if ($_==5 and $test_use_bam==0){&ReverseComplement(); last MYTEST01};
		if ($_==6 and $test_use_bam==0){&Bwa(); last MYTEST01};
#		if ($_==7){&ExtractPaired(); last MYTEST01};
#		if ($_==8){&ExtractUnique(); last MYTEST01};
		if ($_==7){&InsertSizeCollect();last MYTEST01};
		die "Error at $Stage: unknown procedure: $_\n";
	}
}



###Fastq pre-check###################################################
#&PreCheck()
sub PreCheck {
	$cmd='';
	$Stage='Stage '.$stepth.' PreCheck';
	print "\n\n\n###BEGIN###$Stage\n";
	mkdir ("$run_dir", 0766) || die "Error at $Stage: $run_dir\n" unless (-d $run_dir);;
	chdir "$run_dir" || die "Error at $Stage: can not cd dir $run_dir：$!\n";
	my $newdir=$run_dir.'/step'.$stepth.'_fastq';
	mkdir ($newdir, 0766) || die "Error at $Stage: can not create directory $newdir\n";
	chdir "$newdir" or die "Error at $Stage: can not cd dir $newdir：$!\n";
	for (my $i=0; $i<$num_groupofreads;$i++) {
		my $outputR1=$prefix."_L$i"."_R1_precheck.$input_format";
		unlink ($outputR1) if (-e $outputR1);
		my $outputR2=$prefix."_L$i"."_R2_precheck.$input_format";
		unlink ($outputR2) if (-e $outputR2);
		&UnCompress($ReadGroups[$i][0], $outputR1);
		&UnCompress($ReadGroups[$i][1], $outputR2);
		&AddToReadGroup($i, "$newdir/$outputR1", "$newdir/$outputR2", 0);#######                      ###temp
	}
	$stepth++;
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	print "###END###$Stage\n\n\n";
}



###CutAdapt remove short adaptors at 5' or 3' end####################
#&CutAdapt()
sub CutAdapt {
	$cmd='';
	$Stage='Stage '.$stepth.' CutAdapt';
	print "\n\n\n###BEGIN###$Stage\n";
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	my $newdir=$run_dir.'/step'.$stepth.'_cutadapt';
	mkdir ($newdir, 0766) || die "Error at $Stage: can not create directory $newdir\n";
	chdir "$newdir" or die "Error at $Stage: can not cd dir $newdir：$!\n";
	$cmd="$path_cutadapt --format=$input_format --times=$cutadapt_times --overlap=$cutadapt_overlap";
	my %adaptor_list=();
	open (ADAPTORS, "$file_adaptor") || die "Error at $Stage: Can not open the adaptors file $file_adaptor at \n";
	my $adaptor_id='';
	my $adaptor_seq='';
	while (my $line=<ADAPTORS>) {
		chomp $line;
		if ($line=~m/^>(\S+)\s*.*/) {
			if (($adaptor_id ne '') and ($adaptor_seq ne '')){
				if (! exists $adaptor_list{$adaptor_id}) {
					$adaptor_list{$adaptor_id}=$adaptor_seq;
					$adaptor_seq='';
				}
				else {
					die "Error at $Stage: repeated adaptor ID\n";
				}
			}
			$adaptor_id=$1;
		}
		else {
			$adaptor_seq.=$line;
		}
	}
	$adaptor_list{$adaptor_id}=$adaptor_seq;
	close ADAPTORS;
	my @primers_arr=();
	@primers_arr=keys %adaptor_list;
	my ($cmd_r1, $cmd_r2)=($cmd,$cmd);
	my ($test_r1, $test_r2)=(0,0);
	foreach (@primers_arr) {
		if ($_=~/(\S+)_3p$/i) {
			if ($1=~/\S+_r1/) {
				$cmd_r1.=" --adapter $adaptor_list{$_} ";
				$test_r1++;
			}
			elsif ($1=~/\S+_r2/) {
				$cmd_r2.=" --adapter $adaptor_list{$_} ";
				$test_r2++;
			}
			else {
				$cmd.=" --adapter $adaptor_list{$_} ";
				$cmd_r1.=" --adapter $adaptor_list{$_} ";
				$cmd_r2.=" --adapter $adaptor_list{$_} ";
			}
		}
		elsif ($_=~/(\S+)_5p$/i) {
			if ($1=~/\S+_r1/) {
				$cmd_r1.=" --front $adaptor_list{$_} ";
				$test_r1++;
			}
			elsif ($1=~/\S+_r2/) {
				$cmd_r2.=" --front $adaptor_list{$_} ";
				$test_r2++;
			}
			else {
				$cmd.=" --front $adaptor_list{$_} ";
				$cmd_r1.=" --front $adaptor_list{$_} ";
				$cmd_r2.=" --front $adaptor_list{$_} ";
			}
		}
		else {
			$cmd.=" --anywhere $adaptor_list{$_} ";
			$cmd_r1.=" --anywhere $adaptor_list{$_} ";
			$cmd_r2.=" --anywhere $adaptor_list{$_} ";
		}
	}
	for (my $i=0; $i<$num_groupofreads;$i++) {
		my $outputR1=$prefix."_L$i"."_R1_cutadapt.$input_format";
		unlink ($outputR1) if (-e $outputR1);
		my $outputR2=$prefix."_L$i"."_R2_cutadapt.$input_format";
		unlink ($outputR2) if (-e $outputR2);
		if ($test_r1>0) {
			&exec_cmd($cmd_r1." -o $outputR1 $ReadGroups[$i][0]");
		}
		if ($test_r2>0) {
			&exec_cmd($cmd_r2." -o $outputR2 $ReadGroups[$i][1]");
		}
		if ($test_r1==0) {
			&exec_cmd($cmd." -o $outputR1 $ReadGroups[$i][0]");
		}
		if ($test_r2==0) {
			&exec_cmd($cmd." -o $outputR2 $ReadGroups[$i][1]");
		}
		&AddToReadGroup($i, "$newdir/$outputR1", "$newdir/$outputR2", 1);
	}
	$stepth++;
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	print "###END###$Stage\n\n\n";
}



###Trimmomatic trim and clean the reads
#&Trimmomatic()
sub Trimmomatic {
	$cmd='';
	$Stage='Stage '.$stepth.' Trimmomatic';
	print "\n\n\n###BEGIN###$Stage\n";
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	my $newdir=$run_dir.'/step'.$stepth.'_trimmomatic';
	mkdir ($newdir, 0766) || die "Error at $Stage: can not create directory $newdir\n";
	chdir "$newdir" or die "Error at $Stage: can not cd dir $newdir：$!\n";
	$cmd="$path_java -jar $path_trimmomatic PE -threads $num_threads ";
	if (lc($quality_score_format) eq 'phred33') {
		$cmd.=' -phred33 ';
	}
	elsif (lc($quality_score_format) eq 'phred64') {
		$cmd.=' -phred64 ';
	}
	else {
		die "Error at $Stage: Uknown quality score format $quality_score_format\n";
	}
	for (my $i=0; $i<$num_groupofreads;$i++) {
		my $outputR1=$prefix."_L$i"."_R1_trim.$input_format";
		unlink ($outputR1) if (-e $outputR1);
		my $outputR2=$prefix."_L$i"."_R2_trim.$input_format";
		unlink ($outputR2) if (-e $outputR2);
		my $cmdTrimmo=$cmd." $ReadGroups[$i][0] $ReadGroups[$i][1] $outputR1 $prefix"."_L$i".'_R1_trim_unpaired.'."$input_format $outputR2 $prefix"."_L$i".'_R2_trim_unpaired.'."$input_format ILLUMINACLIP:$contaminants:2:30:10 LEADING:18 TRAILING:18 HEADCROP:8 SLIDINGWINDOW:4:18 MINLEN:$trimmo_minlen";
		&exec_cmd($cmdTrimmo);
		&AddToReadGroup($i, "$newdir/$outputR1", "$newdir/$outputR2", 1);
	}	
	$stepth++;
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	print "###END###$Stage\n\n\n";
}



###FastqJoin: join read pairs together if they have $min_overlap overlap for long mate pairs
#&FastqJoin()
sub FastqJoin {
	$cmd='';
	$Stage='Stage '.$stepth.' FastqJoin';
	print "\n\n\n###BEGIN###$Stage\n";
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	my $newdir=$run_dir.'/step'.$stepth.'_fastqjoin';
	mkdir ($newdir, 0766) || die "Error at $Stage: can not create directory $newdir\n";
	chdir "$newdir" or die "Error at $Stage: can not cd dir $newdir：$!\n";
	for (my $i=0; $i<$num_groupofreads;$i++) {
		$cmd="$path_fastqjoin -m $min_overlap ";
		my $outputR1=$prefix."_L$i"."_R1_join.fq";
		unlink ($outputR1) if (-e $outputR1);
		my $outputR2=$prefix."_L$i"."_R2_join.fq";
		unlink ($outputR2) if (-e $outputR2);
		my $output_join=$prefix."_L$i"."_join.fq";
		unlink ($output_join) if (-e $output_join);
		$cmd.=" $ReadGroups[$i][0] $ReadGroups[$i][1] -o $outputR1 -o $outputR2 -o $output_join";
		&exec_cmd($cmd);
		&AddToReadGroup($i, "$newdir/$outputR1", "$newdir/$outputR2", 0);
		if (-s $output_join) {
			$cmd='';
			$cmd="perl -ne \'print length(<>).\"\n\"; <>; <>;\' $output_join > $output_join".'.length';
			&exec_cmd($cmd);
			&sizebin($output_join.'.length', $sizebin, $output_join.'.sizebin');
		}
	}
	$stepth++;
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	print "###END###$Stage\n\n\n";
}



###Trim read to certain max-length
#&FastxTrimmer()
sub FastxTrimmer{
	$cmd='';
	$Stage='Stage '.$stepth.' FastxTrimmer';
	print "\n\n\n###BEGIN###$Stage\n";
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	my $newdir=$run_dir.'/step'.$stepth.'_fastxtrimmer';
	mkdir ($newdir, 0766) || die "Error at $Stage: can not create directory $newdir\n";
	chdir "$newdir" or die "Error at $Stage: can not cd dir $newdir：$!\n";
	$cmd="$path_fastxtrimmer -l $trimmer_minlen";
	for (my $i=0; $i<$num_groupofreads;$i++) {
		my $outputR1=$prefix."_L$i"."_R1_trimmer.$input_format";
		unlink ($outputR1) if (-e $outputR1);
		my $outputR2=$prefix."_L$i"."_R2_trimmer.$input_format";
		unlink ($outputR2) if (-e $outputR2);
		my $cmdR1='';
		$cmdR1=$cmd." -i $ReadGroups[$i][0] -o $outputR1";
		&exec_cmd($cmdR1);
		my $cmdR2='';
		$cmdR2=$cmd." -i $ReadGroups[$i][1] -o $outputR2";
		&exec_cmd($cmdR2);
		&AddToReadGroup($i, "$newdir/$outputR1", "$newdir/$outputR2", 1);
	}
	$stepth++;
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	print "###END###$Stage\n\n\n";
}



###reverse complement
#&ReverseComplement()
sub ReverseComplement {
	$cmd='';
	$Stage='Stage '.$stepth.' ReverseComplement';
	print "\n\n\n###BEGIN###$Stage\n";
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	my $newdir=$run_dir.'/step'.$stepth.'_reversecomplement';
	mkdir ($newdir, 0766) || die "Error at $Stage: can not create directory $newdir\n";
	chdir "$newdir" or die "Error at $Stage: can not cd dir $newdir：$!\n";
	$cmd="$path_fastxreversecomplement";
	if ($quality_score_format eq '-phred33') {
		$cmd.=' -Q33 ';
	}
	for (my $i=0; $i<$num_groupofreads;$i++) {
		my $outputR1=$prefix."_L$i"."_R1_rc.$input_format";
		unlink ($outputR1) if (-e $outputR1);
		my $outputR2=$prefix."_L$i"."_R2_rc.$input_format";
		unlink ($outputR2) if (-e $outputR2);
		my $cmdR1='';
		$cmdR1=$cmd." -i $ReadGroups[$i][0] -o $outputR1";
		&exec_cmd($cmdR1);
		my $cmdR2='';
		$cmdR2=$cmd." -i $ReadGroups[$i][1] -o $outputR2";
		&exec_cmd($cmdR2);
		&AddToReadGroup($i, "$newdir/$outputR1", "$newdir/$outputR2", 1)
	}
	$stepth++;
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	print "###END###$Stage\n\n\n";
}



###BWA mapping#######################################################
#&Bwa()
sub Bwa {
	$cmd='';
	$Stage='Stage '.$stepth.' BWA';
	print "\n\n\n###BEGIN###$Stage\n";
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	if ($test_reIndex) {
		unlink ("ref");
		unlink ("ref.*");
	}
	unless (-l 'ref' or -e 'ref') {
		$cmd="ln -sf $bwaref ref";
		&exec_cmd($cmd);
	};
	unless (-s 'ref.amb' and -s 'ref.ann' and -s 'ref.bwt' and -s 'ref.pac' and -s 'ref.sa') {
		$cmd="$path_bwa index -a is -p ref ref";
		&exec_cmd($cmd);
	}
	my $newdir=$run_dir.'/step'.$stepth.'_bwa';
	mkdir ($newdir, 0766) || die "Error at $Stage: can not create directory $newdir\n";
	chdir "$newdir" or die "Error at $Stage: can not cd dir $newdir：$!\n";
	for (my $i=0; $i<$num_groupofreads;$i++) {
		$cmd="$path_bwa aln -t $num_threads $run_dir".'/ref ';
		my $outputR1=$prefix."_L$i"."_R1_bwaaln.sai";
		unlink ($outputR1) if (-e $outputR1);
		my $outputR2=$prefix."_L$i"."_R2_bwaaln.sai";
		unlink ($outputR2) if (-e $outputR2);
		my $output_bam=$prefix."_L$i"."_bwaalnQ$min_mapq.bam";
		unlink ($output_bam) if (-e $output_bam);
		my $output_sortbam=$prefix."_L$i"."_R2_bwaalnQ$min_mapq.sort.bam";
		my $output_sortbam_prefix=$prefix."_L$i"."_bwaalnQ$min_mapq.sort";
		unlink ($output_sortbam) if (-e $output_sortbam);
		my $cmdR1='';
		$cmdR1=$cmd." $ReadGroups[$i][0] > $outputR1";
		&exec_cmd($cmdR1);
		my $cmdR2='';
		$cmdR2=$cmd." $ReadGroups[$i][1] > $outputR2";
		&exec_cmd($cmdR2);
		if (-s $outputR1 and -s $outputR2) {
#			my $cmd_sampe="$path_bwa sampe -a $maximum_insertsize -r \'".'@RG\tID:L'.$i.'\tSM:L'.$i."\'".' '.$run_dir.'/ref '."$outputR1 $outputR2 $ReadGroups[$i][0] $ReadGroups[$i][1]".' | '.$path_samtools.' view -bS -q '.$min_mapq.' -F 12 - > '.$output_bam;
			my $cmd_sampe="$path_bwa sampe -a $maximum_insertsize ".$run_dir.'/ref '."$outputR1 $outputR2 $ReadGroups[$i][0] $ReadGroups[$i][1]".' | '.$path_samtools.' view -bS -q '.$min_mapq.' -F 12 - > '.$output_bam;
			&exec_cmd($cmd_sampe);
			$cmd='';
			$cmd="$path_samtools sort $output_bam $output_sortbam_prefix";
			&exec_cmd($cmd);
			$cmd='';
			$cmd="$path_java -Xmx$picard_memory -jar $path_picard/MarkDuplicates.jar I=$output_sortbam O=$output_sortbam_prefix.rmdup.bam REMOVE_DUPLICATES=TRUE METRICS_FILE=$output_sortbam_prefix.rmdup.log ASSUME_SORTED=TRUE";
			&exec_cmd($cmd);
		}
		else {
			die "Error at $Stage: BWA aln output error\n";
		}
		if (-s "$output_sortbam_prefix.rmdup.bam") {
			$BamGroups[$i]=$newdir.'/'.$output_sortbam_prefix.'.rmdup.bam';
		}
		else {
			die "Error at $Stage: BWA sampe output error\n";
		}
	}
	$stepth++;
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	print "###END###$Stage\n\n\n";
}



sub InsertSizeCollect {
	$cmd='';
	$Stage='Stage '.$stepth.' InsertSizeCollect';
	print "\n\n\n###BEGIN###$Stage\n";
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	my $newdir=$run_dir.'/step'.$stepth.'_InsertSizeCollect';
	mkdir ($newdir, 0766) || die "Error at $Stage: can not create directory $newdir\n";
	chdir "$newdir" or die "Error at $Stage: can not cd dir $newdir：$!\n";
	for (my $i=0; $i<$num_groupofreads;$i++) {
		$cmd="perl $Bin/utils/getInsertSize_from_sambam_LUF.pl --min_insertsize $minimum_insertsize --max_insertsize $maximum_insertsize --min_mapq $min_mapq --path_samtools $path_samtools --input $BamGroups[$i] ";
		my $outputPaired=$prefix."_L$i"."_Paired.sum";
		unlink ($outputPaired) if (-e $outputPaired);
		my $outputUnique=$prefix."_L$i"."_Unique.sum";
		unlink ($outputUnique) if (-e $outputUnique);
		my $cmdPaired=$cmd." --output $outputPaired ";
		&exec_cmd($cmdPaired);
		my $cmdUnique=$cmd." --output $outputUnique --tag XT:A:U";
		&exec_cmd($cmdUnique);
	}
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	print "###END###$Stage\n\n\n";
}



=ignore
###Extract read pairs from @BamGroups and re-map to reference
#&ExtractPaired()
sub ExtractPaired {
	$cmd='';
	$Stage='Stage '.$stepth.' ExtractPaired';
	print "\n\n\n###BEGIN###$Stage\n";
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	my $newdir=$run_dir.'/'.$stepth.'.extractpaired';
	mkdir ($newdir, 0766) || die "Error at $Stage: can not create directory $newdir\n";
	chdir "$newdir" or die "Error at $Stage: can not cd dir $newdir：$!\n";
	$cmd="$path_java' -Xmx'$picard_memory' -jar '$path_picard'/SamToFastq.jar '";
	for (my $i=0; $i<$num_groupofreads;$i++) {
		my $outputR1=$prefix."_L$i"."_paired_R1.$input_format";
		unlink ($outputR1) if (-e $outputR1);
		my $outputR2=$prefix."_L$i"."_paired_R2.$input_format";
		unlink ($outputR2) if (-e $outputR2);
		$cmd.=" 'INPUT='$BamGroups[$i]' FASTQ='$outputR1' SECOND_END_FASTQ='$outputR2' INCLUDE_NON_PRIMARY_ALIGNMENTS=false'";
		&exec_cmd($cmd);
		&AddToReadGroup($i, "$newdir/$outputR1", "$newdir/$outputR2", 1);
	}
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	print "###END###$Stage\n\n\n";
	&Bwa();
}



###Extract Unique mapped read pairs from @BamGroups and re-map to reference
#&ExtractUnique()
sub ExtractUnique {
	$cmd='';
	$Stage='Stage '.$stepth.' ExtractUnique';
	print "\n\n\n###BEGIN###$Stage\n";
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	my $newdir=$run_dir.'/'.$stepth.'.extractunique';
	mkdir ($newdir, 0766) || die "Error at $Stage: can not create directory $newdir\n";
	chdir "$newdir" or die "Error at $Stage: can not cd dir $newdir：$!\n";
	$cmd="$path_java' -Xmx'$picard_memory' -jar '$path_picard'/SamToFastq.jar '";
	for (my $i=0; $i<$num_groupofreads;$i++) {
		my $outputR1=$prefix."_L$i"."_unique_R1.$input_format";
		unlink ($outputR1) if (-e $outputR1);
		my $outputR2=$prefix."_L$i"."_unique_R2.$input_format";
		unlink ($outputR2) if (-e $outputR2);
		my $output_filterBam=$prefix.'_L'.$i.'_R2_bwaalnQ'.$min_mapq.'.sort.rmdup.bam';
		unlink ($output_filterBam) if (-e $output_filterBam);
		my $cmd_unique="$path_samtools' view -h '$BamGroups[$i]' | egrep '\"(^@|XT\:A\:U)\"' | '$path_samtools' view -bS - > '$output_filterBam";
		&exec_cmd($cmd_unique);
		$cmd.=" 'INPUT='$output_filterBam' FASTQ='$outputR1' SECOND_END_FASTQ='$outputR2' INCLUDE_NON_PRIMARY_ALIGNMENTS=false'";
		&exec_cmd($cmd);
		&AddToReadGroup($i, "$newdir/$outputR1", "$newdir/$outputR2", 1);
	}
	chdir "$run_dir" or die "Error at $Stage: can not cd dir $run_dir：$!\n";
	print "###END###$Stage\n\n\n";
	&Bwa();
}
=cut



###Add path to a file
#&AddFilePath(files)
sub AddFilePath {
	my @AFPfiles=@_;
	use Cwd;
	my $AFPcurrentdir=getcwd;
	for (my $AFPnum_file=0;$AFPnum_file<scalar(@AFPfiles);$AFPnum_file++) {
		unless ($AFPfiles[$AFPnum_file]=~/.*\//) {
			$AFPfiles[$AFPnum_file]="$AFPcurrentdir/$AFPfiles[$AFPnum_file]";
		}
	}
	return (@AFPfiles);
}



### Backup file
#&backup(file)
sub backup {
	my $BUori_name=shift @_;
	my $BU_new_name='';
	$BU_new_name=$BUori_name.".bak.".int(rand(10000));
	if (-e "$BU_new_name") {
		&backup($BUori_name);
	}
	rename($BUori_name, $BU_new_name);
	print "The Original file $BUori_name was renamed to $BU_new_name\n";
}



###determine how to process those files existed
#&OutputExistsTest(file)
sub OutputExistsTest {
	my $FET_fileid=shift @_;
	chomp $FET_fileid;
	if (-e $FET_fileid) {
		if ($verbose) {
			print "The $FET_fileid file specified already existed\nNeed to overwrite [y] or backup [n] or others [AnyKey], default[y]:\n";
			my $FET_test=''; $FET_test=<STDIN>;
			chomp $FET_test;
			if (lc ($FET_test) eq ('y'||'yes')) {
				unlink($FET_fileid);
			}
			elsif (lc ($FET_test) eq ('n'||'no')) {
				&backup($FET_fileid);
			}
			else {
				die "Please specify a new name for output\n";
			}
		}
		else {
			unlink($FET_fileid);
		}
	}
}



###Tect if a command exists
#&TestCmdExist(path/to/cmd)
sub TestCmdExist {
	my $TCEfile=shift @_;
	if (-f $TCEfile) {
		if (-e $TCEfile and -x _) {
			return 1;
		}
		else {
			die "File $TCEfile exists, but not executable\n";
		}
	}
	elsif (-d $TCEfile) {
		die "$TCEfile is not a file name, but a directory name\n";
	}
	else {
		die "Can not find $TCEfile\n";
	}
}



###Return: YYYYMMDD hh:mm:ss
#&mytime()
sub mytime {
	my($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime();		$year += 1900;
	$mon  += 1;
	my $btime = sprintf("%04d%02d%02d %02d:%02d:%02d",$year,$mon,$day,$hour,$min,$sec);
	return $btime;
}



###Process command
#&exec_cmd(cmd)
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
#		print "Error, cmd: $cmd died with ReturnCode $return_code\n";
		die "Error, cmd: $cmd died with ReturnCode $return_code\n";
		return $return_code;
	}
	else {
		print "Finished command: $cmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n";
		return $return_code;
	}
}



###Put sizes in a bin
#&sizebin(INPUTFILE, $bin, OUTPUTFILE);
sub sizebin {
	use Scalar::Util qw(looks_like_number);
	my ($input, $bin, $output)=@_;
	my %count=();
	my $min=0;
	my $max=0;
	open (INPUT, "$input") || die "Can not open input file\n";
	while (my $line=<INPUT>) {
		next if ($line=~/^#|^\*/);
		chomp $line;
		my @arr=();
		@arr=split(/\s+|,|;/,$line);
		my $col=0;
		for (my $i=0; $i<scalar(@arr);$i++){
			if (looks_like_number($arr[$i])) {
				$col=$i;
				last;
			}
		}
		if (scalar(@arr)>1 and $col < scalar(@arr)-1 and looks_like_number($arr[$col+1])) {
			$count{$bin*int($arr[$col]/$bin)}+=$arr[$col+1];
		}
		else {
			$count{$bin*int($arr[$col]/$bin)}++;
		}
		$min=$arr[$col] if ($arr[$col]<$min);
		$max=$arr[$col] if ($arr[$col]>$max);
	}
	close INPUT;
	my @countkey=keys %count;
	@countkey=sort {$a <=> $b} @countkey;
	unlink ($output) if (-e $output);
	open (OUTPUT, ">>$output") || die "Error at $Stage: can not output to file $output";
	print OUTPUT "Min: $min\tMax:$max\tBin:$bin\n";
	for (my $j = $countkey[0]; $j<$max; $j+=$bin) {
		if (exists $count{$j}) {
			print OUTPUT $j."\t".$count{$j}."\n";
		}
		else {
			print OUTPUT $j."\t0\n";
		}
	}
	close OUTPUT;
	no Scalar::Util;
}



###Push Read to @ReadGroups
#&AddToReadGroup($i, $R1, $R2, run_fastqc)
sub AddToReadGroup {
	my ($i, $R1, $R2, $run_this_fastqc)=@_;
	if (-s $R1 and -s $R2) {
		@{$ReadGroups[$i]}=($R1, $R2);
		if ($run_this_fastqc){
			&FastQC($R1);
			&FastQC($R2);
		}
	}
	else {
		die "Error at $Stage: output empty/none files\n";
	}
}



###Push bam to @BamGroups
#&AddToBamGroup($i, $bam)
sub AddToBamGroup {
	my ($i, $bam)=@_;
	if (-s $bam) {
		$BamGroups[$i]=$bam;
	}
	else {
		die "Error at $Stage: output empty/none files\n";
	}
}



###produce fastqc report file (.zip) 
#&FastQC(Fastq_File)
sub FastQC {
	my $fastq=shift @_;
	&exec_cmd("$path_fastqc $fastq -o . -f fastq -t $num_threads --noextract") if ($test_run_fastqc>0);
}



###UnCompress
#UnCompress (Compress_file_In, UnCompressedFile_Out)
sub UnCompress {
	my ($in, $out)=@_;
	print "UnCompress $in to $out ...\n";
	MYTEST: {
		if ($in=~/\.fastq$|\.fq$/i) {&exec_cmd("ln -sf $in $out"); last MYTEST;}
		if ($in=~/\.tar\.gz$/i) {&exec_cmd("tar -xvf --to-stdout $in > $out"); last MYTEST;}
		if ($in=~/\.gz$|\.gzip$/i) {&exec_cmd("gunzip -c $in > $out"); last MYTEST;}
		if ($in=~/\.bz2$|\.bzip2$/i) {&exec_cmd("tar -xjf --to-stdout $in > $out"); last MYTEST;}
	}
}



###Get environment variables
#&GetEnvPath(VariableName)
sub GetEnvPath {
	my $varable=shift @_;
		unless (defined $path_picard) {
		if ($ENV{'PICARD_HOME'}) {
			$path_picard=$ENV{'PICARD_HOME'};
		}
		else {
			die "Error at $Stage: can not find PICARD_HOME or path to PICARD\n";
		}
	}
}
