#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use constant USAGE =><<EOH;

SYNOPSIS:

perl $0 [Options]
Version: LUFUHAO20140801

Descriptions:
	Calculate insert size for paired end or mate paired NGS. 
	Input is sorted bam.
	Require: samtools

Options:
	--help/-h
		Print this help/usage;
	--input|-i <File>
		[Msg] Sorted bam file;
	--output|-o <File>
		[Msg] Output file;
	--min_insertsize <INT>
		[Opt] Minimum insert size threshold to keep; Default: 0;
	--max_insertsize <INT>
		[Opt] Maximum insert size threshold to keep; Default: unset;
	--path_samtools </Path/to/samtools>
		[Opt] path/to/samtools if it's not in your PATH;
	--min_mapq <INT>
		[Opt] Minimum mapping quality threshold to keep an alignment;
	--tag|-g <STR>
		[Opt] Sam flags used to keep an alignment;
	--verbose
		[Opt] Detailed output for trouble-shooting;
	--version/-v
		[Opt] Print current SCRIPT version;

Example:
	perl $0 \
	--input my.sorted.bam \
	-o insert.size.txt \
	--min_insertsize 20000 \
	--max_insertsize 60000 \
	--path_samtools /usr/local/samtools/v0.1.18/bin/samtools \
	--min_mapq 5 \
	--tag XT:A:U

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
our ($help, $input, $verbose, $version);
our ($output, $min_insertsize, $max_insertsize, $path_samtools, $min_mapq, $tag);
GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"output|o:s" => \$output,
	"min_insertsize:i" => \$min_insertsize,
	"max_insertsize:i" => \$max_insertsize,
	"path_samtools:s" => \$path_samtools,
	"min_mapq:i" => \$min_mapq,
	"tag|g" => \$tag,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
#die USAGE unless (@ARGV);
($help or $version) and die USAGE;



###Configure and default#############################################
$path_samtools='samtools' unless (defined $path_samtools);
$min_mapq=0 unless (defined $min_mapq);
$min_insertsize=0 unless (defined $min_insertsize);
###Input and output##################################################
our $cmd='';
my $file_basename='';
if ($input=~/(.*)\.sam$/i) {
	$file_basename=$1;
	$cmd="$path_samtools view -bS -f 3 -F 12 $input > $file_basename.bam";
	&exec_cmd($cmd);
	$cmd="$path_samtools sort $file_basename.bam $file_basename.sort";
	&exec_cmd($cmd);
	unlink("$file_basename.bam");
	open (INPUT, "samtools view -h $file_basename.sort.bam |") || die "Can not view this file using samtools\n";
	print "Sam format INPUT accepted: $input\n" if (defined $verbose);
}
elsif ($input=~/(.*)\.bam$/i) {
	$file_basename=$1;
	unlink ("$input.bai") if (-e "$input.bai");
	if (&exec_cmd_return("$path_samtools index $input")) {
		$cmd="$path_samtools sort $input $file_basename.sort";
		&exec_cmd($cmd);
		open (INPUT, "samtools view -h $file_basename.sort.bam |") || die "Can not view this file using samtools\n";
		print "Unsorted Bam format INPUT accepted: $input\n" if (defined $verbose);
	}
	else {
		open (INPUT, "samtools view -h $input |") || die "Can not view this file using samtools\n";
		print "Sorted Bam format INPUT accepted: $input\n" if (defined $verbose);
	}
}
&OutputExistsTest($output);
our @tags=split(/,/, $tag) if (defined $tag);



###Main program######################################################
our %insertsize=();
our %id=();
our @ref=();
my $temp01=0;
my $ref='';
my $previous_ref='';
while (my $line=<INPUT>) {
	chomp $line;
	if ($line=~m/^\@/) {
		if ($line=~m/^\@SQ\s+SN:(\S+)\s+LN:\d+/) {
			push (@ref, $1);
		}
		next;
	}
	my @arr=split(/\t/, $line);
#	print "$arr[2]\n";                             ###test###
	die "Wrong Bam input $line\n" if (@arr < 11);
#0	1	2	3	4	5	6	7	8	9	10
#read	FLAG	ref	pos	mapq	cigar	mate	pos	Tlen	seq	qual
#3NG5HQ1:156:C2C0AACXX:8:1302:8229:56841 99      1018_2al        1       27      16S68M7S        =       106     212     TGTGAGACTATCTTCTCCTTTTTGTCTTCTCCACAACCACCATTCTATTCCACCTATAGTGCTATATCCATGGCTCACGCTCATGTATTGC     HHHJGIIJJJJJJJJJJJJJJJJGHGIGEFHIHGGGIIIIGJJJIJJIJJJIJIJJJIIHIJJIJEHIHHHHHHFFFFFEDDDDDDDEDD>     RG:Z:ParRoot1   MD:Z:68 NH:i:2  HI:i:1  NM:i:0  SM:i:27 XQ:i:40 X2:i:34 XO:Z:CM PG:Z:M
	next if ((! $arr[1] & 0x3) or ($arr[1] & 0xc) or (! defined $arr[2]) or $arr[3]<0 or $arr[4]<$min_mapq or $arr[8]==0);#($arr[6] ne '*' or $arr[6] ne '=') or $arr[7]>0
	if (defined $tag) {
		foreach my $ind_tag (@tags) {
			next unless ($line=~m/$ind_tag/);
		}
	}
	(my $base_name =$arr[0])=~s/\/.*//;
#	print "$arr[2]\n";                             ###test###
#	$base_name=quotemeta($base_name);
	$ref=$arr[2];
	if ($temp01==0) {
		$previous_ref=$ref;
		$temp01=1;
	}
	if ($ref eq $previous_ref) {
		if (defined $id{$base_name}) {
			if (defined $insertsize{$base_name}) {
				$insertsize{$base_name}='Multiple';
			}
			elsif (scalar (@{$id{$base_name}})==1) {
				@{$id{$base_name}[1]}=($arr[2], $arr[3], length($arr[9]));
			}
			if (scalar (@{$id{$base_name}})==2) {
				&CalculateInsertSize($base_name, @{$id{$base_name}[0]}, @{$id{$base_name}[1]});
				delete $id{$base_name};
			}
		}
		else {
			@{$id{$base_name}[0]}=($arr[2], $arr[3], length($arr[9]));
		}
	}
	else {
		&WriteInsertSize($previous_ref);
		@{$id{$base_name}[0]}=($arr[2], $arr[3], length($arr[9]));
		$previous_ref=$ref;
	}
}
&WriteInsertSize($ref);
close INPUT;



#####################################################################
### sub function ####################################################
#####################################################################

###Current time
#&mytime()
sub mytime {
	my($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime();	$year += 1900;
	$mon  += 1;
	my $btime = sprintf("%04d%02d%02d %02d:%02d:%02d",$year,$mon,$day,$hour,$min,$sec);
	return $btime;}



###Process command
#&exec_cmd(cmd)
sub exec_cmd {
	my ($ECcmd) = @_;
	print &mytime()."CMD: $ECcmd\n"  if (defined $verbose);
	my $start_time = time();
	my $return_code = system($ECcmd);
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
#		print "Error, cmd: $ECcmd died with ReturnCode $return_code\n";
		die "Error, cmd: $ECcmd died with ReturnCode $return_code\n";
		return $return_code;
	}
	else {
		print "Finished command: $ECcmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n" if (defined $verbose);
		return $return_code;
	}
}



sub exec_cmd_return {
	my ($RCRcmd) = @_;
	print &mytime()."CMD: $RCRcmd\n" if (defined $verbose);
	my $start_time = time();
	my $return_code = system($RCRcmd);
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
		print "Error, cmd: $RCRcmd died with ReturnCode $return_code\n" if (defined $verbose);
#		die "Error, cmd: $RCRcmd died with ReturnCode $return_code\n";
		return $return_code;
	}
	else {
		print "Finished command: $RCRcmd\nat ".&mytime()."\nRunning time:(".($end_time - $start_time)." seconds) with Returncode: $return_code\n" if (defined $verbose);
		return $return_code;
	}
}



###Calculate insert size
#&CalculateInsertSize(ref1, pos1, len1, ref2, pos2, len2)
sub CalculateInsertSize {
	my ($read_id, $ref1, $pos1, $len1, $ref2, $pos2, $len2)=@_;
	my $CISinsertsize=0;
#	print "$read_id, $ref1, $pos1, $len1, $ref2, $pos2, $len2\n";          ###for test###
#	$read_id=quotemeta($read_id);
#	print "$pos1\t$pos2\n$len2";                                           ###for test###
	$CISinsertsize=$pos2-$pos1;
#	print "$CISinsertsize\n";                                              ###for test###
	$CISinsertsize=abs($CISinsertsize)+$len2+1;
#	print "$CISinsertsize\n";                                              ###for test###
	if (defined $max_insertsize) {
		$insertsize{$read_id}=$CISinsertsize if (($CISinsertsize>$min_insertsize) and ($CISinsertsize<$max_insertsize));
	}
	else {
		$insertsize{$read_id}=$CISinsertsize if ($CISinsertsize>$min_insertsize);
	}
#	print $insertsize{$read_id}."\n";                                      ###for test###
}



###WriteInsertSize
#&WriteInsertSize(ref)
sub WriteInsertSize {
	my $ref_id=shift @_;
	my @read_ids=keys %insertsize;
#	print scalar(@read_ids)."\n";                                          ###for test###
#	print "$ref_id\n";                                                     ###for test###
	if (defined $output) {
		open (OUTPUT, ">>$output") || die "Can not write into file: $output\n";
		print OUTPUT "# $ref_id\n";
		foreach my $id (@read_ids) {
			if (defined $id and ($id ne '')) {
				print OUTPUT $id."\t".$insertsize{$id}."\n";
			}
		}
		close OUTPUT;
	}
	else {
		print "# $ref_id\n";
		foreach my $id (@read_ids) {
			if (defined $id and ($id ne '')) {
				print $id."\t".$insertsize{$id}."\n";
			}
		}
		close OUTPUT;
	}
	%insertsize=();
	%id=();
}



###determine how to process those files existed
#OutputExistsTest($file)
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
			&backup($FET_fileid);
		}
	}
}
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
