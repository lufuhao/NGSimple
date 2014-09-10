#!/bin/sh
#Require zlib, GNU ncurses library
#APT
#apt-get install zlib1g zlib1g-dev
#apt-get install libncurses5 libncurses5-dev
#YUM
#yum install zlib zlib-devel
#yum install ncurses ncurses-devel
RunDir=$(cd `dirname $0`; pwd)
MachType=$(uname -m)

if [ ! -d $RunDir/samtools ]; then
mkdir -p $RunDir/samtools
fi
cd $RunDir/samtools
if [ ! -s $RunDir/samtools/samtools*.tar.bz2 ]; then
  echo "Please download samtools source package to $RunDir/samtools/ from website: http://sourceforge.net/projects/samtools"
  exit 0
fi
samtools_package=`ls samtools*.tar.bz2`
samtools_basename=${samtools_package%.*.*}
samtools_version=${samtools_basename#*-}
if [ -d $RunDir/samtools/v$samtools_version ]; then
  rm -rf $RunDir/samtools/v$samtools_version
fi
mkdir -p $RunDir/samtools/v$samtools_version/src
mv $samtools_package $RunDir/samtools/v$samtools_version/src/
cd $RunDir/samtools/v$samtools_version/
tar xvf $RunDir/samtools/v$samtools_version/src/$samtools_package > /dev/null
mv samtools* $MachType
cd $RunDir/samtools/v$samtools_version/$MachType
make CXXFLAGS+=" -fPIC" > make.log
if [ ! -s ./samtools ]; then
  echo "Compiling Samtools failed"
  exit 0
fi
mkdir -p 'bin' 'include' 'lib' 'man/man1'
cp samtools misc/ace2sam misc/maq2sam-long misc/maq2sam-short misc/md5fa misc/md5sum-lite misc/wgsim misc/blast2sam.pl misc/bowtie2sam.pl misc/export2sam.pl misc/interpolate_sam.pl misc/novo2sam.pl misc/plot-bamstats misc/psl2sam.pl misc/sam2vcf.pl misc/samtools.pl misc/seq_cache_populate.pl misc/soap2sam.pl misc/varfilter.py misc/wgsim_eval.pl misc/zoom2sam.pl $RunDir/samtools/v$samtools_version/$MachType/bin/
cp ./bcftools/bcftools $RunDir/samtools/v$samtools_version/$MachType/bin/
cp ./bcftools/vcfutils.pl $RunDir/samtools/v$samtools_version/$MachType/bin/
cp *.h $RunDir/samtools/v$samtools_version/$MachType/include/
cp ./bcftools/*.h $RunDir/samtools/v$samtools_version/$MachType/include/
cp *.a $RunDir/samtools/v$samtools_version/$MachType/lib/
cp ./bcftools/libbcf.a $RunDir/samtools/v$samtools_version/$MachType/lib/
cp samtools.1 $RunDir/samtools/v$samtools_version/$MachType/man/man1/
if [ -d $RunDir/bin ]; then
  mkdir -p $RunDir/bin
fi
ln -sf $RunDir/samtools/v$samtools_version/$MachType/bin/samtools $RunDir/bin/
echo "*******************************************"
echo "************ Important ********************"
echo "*******************************************"
echo "Add to /etc/profile or ~/.bashrc\n\n"
echo "###SAMtools"
echo "export PATH=$RunDir/samtools/v$samtools_version/$MachType/bin:$RunDir/samtools/v$samtools_version/$MachType/misc"':$PATH'
echo "export CPLUS_INCLUDE_PATH=$RunDir/samtools/v$samtools_version/$MachType/include"':$CPLUS_INCLUDE_PATH'
echo "export C_INCLUDE_PATH=$RunDir/samtools/v$samtools_version/$MachType/include"':$C_INCLUDE_PATH'
echo "export LIBRARY_PATH=$RunDir/samtools/v$samtools_version/$MachType/lib"':$LIBRARY_PATH'
echo "export MAN_PATH=$RunDir/samtools/v$samtools_version/$MachType/man"':$MAN_PATH'
echo "\n\n\n"

exit 1
