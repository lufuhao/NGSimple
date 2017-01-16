#NGSimple (Next-Generation illumina SIMulation PipeLinE) 
------------------------------------------------------------------
## 1. Usage

This pipeline, mainly written in perl is designed to simulate NGS 
parameters (illumina sequencing type, read length, insert size and
 coverage) based on a given reference, and give a report to choose
 the best parameter combination that could contribute to a good 
assembly.

**simulation.pl** generates simulated reads from a reference in fasta
 format. It's recommended to *start with a short reference* (~5M 
bases), in case that the number of generated reads is too large for
 assembly.

    input: reference + parameter

    output: ./1.fastq/reads.fastq

**assemly.pl** assemble the reads (./1.fastq/reads.fastq) from those 
designed libraries to assemblies. based on the parameters specified 
in *CONFIG* file. It will generate assmblies with different K-mer 
(./2.assembly/*), and choose the assembly with largest n50. Finally,
 it will produce a quartile summary and mummerplot for each for you 
to select a good parameter combination.

**insert_size.pl** will take real Illumina data as input, trim/clean (
trimmomatic, evalaute the proportion of overlapped pairs (fastq-join),
 reverse complement (fastx_toolkit) , align to reference sequence 
(bwa) and then determine the insert size distribution. You could 
design the steps to ingnore any step in case it's not quite necessary. 
In addition, you could provide any SAM flags to keep an alignment, 
for example, if read mapping is produced by BWA, XT:A:U flag could be
 used to filter the unique mapped reads.

NOTE: the insert size is NOT collected by PICARD CollectInsertSizeMetrics.jar

----------------------------------------------------------------------

##2. Installation

This pipeline makes use of several external programs. So you may need to 
separately install them before running this pipeline. The bash scripts in
 ProgramRoot/utils/ might help you install them if you did not install or 
they are not available in PATH. 

###2.1 Some necessary linux dependancies

*    make, cmake, gcc, g++, gcc-c++, tar, unzip, GD, etc
*    zlib, boost, ncurses, sparsehash, etc
*    PerlIO::gzip, GD::Graphics::bars, Getopt::Long, List::Util, Cwd, FindBin, etc

###2.2 External programs:

**Program** | **Necessary/Optional**
----------- | ----------------------
[Mason](https://www.seqan.de/projects/mason/) | Msg
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | Msg
[Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) | Msg
[MUMmer](http://mummer.sourceforge.net/) | Msg
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | Msg
[CutAdapt](https://code.google.com/p/cutadapt/) | Opt
[fastq-join](https://code.google.com/p/ea-utils/) | Opt
[fastx_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html) | Opt
[BWA](http://bio-bwa.sourceforge.net/) | Msg
[SAMtools](http://samtools.sourceforge.net/) | Msg
[PICARD_tools](http://picard.sourceforge.net/) | Msg

NOTE: the programs need to set the environmental variables (PATH, 
CPLUS_INCLUDE_PATH, LD_LIBRARY_PATH, etc) to be easily found. 

The pipeline ifself does NOT need compilation. Just uncompress it somewhere
 you want, execuate the /path/to/program_root/ngsimple. Put in PATH if 
you want (/etc/profile or ~/.bashrc): 
	export PATH=/path/to/program_root:$PATH


##3. Running

###3.1 configuring parameters

CONFIG file

	################## File border ######################################
	#Example:
	#reference=/Full/path/to/reference.fasta
	reference=
	#example:
	#tab-delimited
	#LibID	libtype	readlength	libcov	insert_size
	#lib1	pe	200	30	600
	#lib2	se	300	5	0
	#lib3	mp	150	10	7000
	lib1
	lib2
	lib3
	#Example:
	#UniqueFolderName: assembly name (space free)
	#UniqueFolderName	index combination
	#tab-delimited
	*IMPORTANT: should support THREE libraries each folder at most;
	#Group1	lib1	lib2
	*will assemble lib1 and lib2 into one assembly under the folder Group1
	#Group2	lib1	lib2	lib3
	*will assemble lib1 and lib2 into one assembly under the folder Group1
	Assem1
	Assem2
	Assem3
	################## File border ######################################


###3.2 Running NGSimple

bash$ ngsimple -i CONFIG -t 5

NOTE: Since Velvet is extremely memory-consuming, it is highly recommended 
that users estimate how much memory will be needed by velvetg before 
submitting their jobs. The estimation may be calcuated using the following 
formula:

Required RAM (MB) = (-109635 + 18977*ReadSize + 86326*GenomeSize + 233353*NumReads - 51092*K)/1024

+ Read size is in bases.

+ Genome size is in millions of bases (Mb)

+ Number of reads is in millions

+ K is the kmer hash value used in velveth

From: http://seqanswers.com/forums/showthread.php?t=2101


###3.3 Determine insert size distribution (NOT using Picard CollectInsertSizeMetrics.jar)
	bash$ 3.insertsize.pl --help

##4. Trouble-shooting
Every perl script could be separately used for advanced users. In case
 you have any problem or question about this pipeline, please E-mail 
me at <Fu-Hao.Lu@jic.ac.uk>


##Author:
  Fu-Hao Lu
  Post-Doctoral Scientist in Micheal Bevan laboratory
  Cell and Developmental Department, John Innes Centre
  Norwich NR4 7UH, United Kingdom
  E-mail: <Fu-Hao.Lu@jic.ac.uk>
