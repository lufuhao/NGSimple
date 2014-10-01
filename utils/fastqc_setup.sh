#!/bin/sh
#Require: JAVA (>1.6.0), unzip
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)

if which fastqc 2>/dev/null; then
  echo "FastQC already installed"
  exit 0
fi

if [ ! -d $RunDir/fastqc ]; then
  mkdir -p $RunDir/fastqc
fi
cd $RunDir/fastqc
if [ ! -s fastqc_*.zip ]; then
  echo "Trying to download FASTQC to $RunDir/fastqc/ from website: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip"
  wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip
fi
if [ ! -s fastqc_*.zip ]; then
  echo "FastQC package is not present and can not bed downloaded"
  exit 1
fi
fastqc_package=`ls fastqc_*.zip`
fastqc_basename=${fastqc_package%.zip}
fastqc_version=${fastqc_basename#fastqc_}
if [ "v$fastqc_version" = "v" ]; then
  echo "Can not detect FastQC version"
  exit 1
fi
echo "FASTQC version $fastqc_version detected"
if [ -d $RunDir/fastqc/$fastqc_version ]; then
  rm -rf $RunDir/fastqc/$fastqc_version
fi
mkdir -p $RunDir/fastqc/$fastqc_version/src
mv $fastqc_package $RunDir/fastqc/$fastqc_version/src/
cd $RunDir/fastqc/$fastqc_version/
unzip $RunDir/fastqc/$fastqc_version/src/$fastqc_package >/dev/null
if [ ! -d FastQC ]; then
  echo "Failed to un-compress $RunDir/fastqc/$fastqc_version/src/$fastqc_package"
  exit 1
fi
mv FastQC $MachType
cd $RunDir/fastqc/$fastqc_version/$MachType
echo "Installing FastQC ..."
if [ ! -s fastqc ]; then
  echo "unzip FastQC failed"
  exit 1
fi
chmod 755 ./fastqc
if [ -d $RunDir/bin ]; then
  mkdir -p $RunDir/bin
fi
ln -sf $RunDir/fastqc/$fastqc_version/$MachType/fastqc $RunDir/bin/
echo "*******************************************"
echo "************ Important ********************"
echo "*******************************************"
echo "Add to /etc/profile or ~/.bashrc\n\n"
echo "###FastQC"
echo "export PATH=$RunDir/fastqc/$fastqc_version/$MachType"':$PATH'
echo "\n\n\n"

exit 0
