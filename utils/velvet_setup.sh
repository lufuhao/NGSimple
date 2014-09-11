#!/bin/sh
#Require: zlib, gcc, (64bit), openmpi
###APT
#apt-get install libopenmpi1.6 libopenmpi-dev openmpi-bin openmpi-common
###YUM
#yum install openmpi openmpi-devel
RunDir=$(cd `dirname $0`; pwd)
MachType=$(uname -m)

if which velvetg 2>/dev/null; then
  echo "velvet already installed"
  exit 0
fi

if [ -d $RunDir/velvet ]; then
  mkdir -p $RunDir/velvet
fi
cd $RunDir/velvet
if [ -d $RunDir/velvet/version ]; then
  rm -rf $RunDir/velvet/version
fi
mkdir -p $RunDir/velvet/version/src
if [ -s velvet_*.tgz ]; then
  echo "velvet package "`ls velvet_*.tgz`" is detected"
  mv velvet_*.tgz $RunDir/velvet/version/src/
  cd $RunDir/velvet/version/
  tar xvf $RunDir/velvet/version/src/velvet_*.tgz > /dev/null
  mv velvet_* velvet
else
  echo "Trying to download latest velvet to $RunDir/velvet/version/src/ from http://www.ebi.ac.uk/~zerbino/velvet/velvet_latest.tgz"
  wget http://www.ebi.ac.uk/~zerbino/velvet/velvet_latest.tgz
  if [ -s velvet_*.tgz ]; then
    echo "velvet package "`ls velvet_*.tgz`" is detected"
    mv velvet_*.tgz $RunDir/velvet/version/src/
    cd $RunDir/velvet/version/
    tar xvf $RunDir/velvet/version/src/velvet_*.tgz > /dev/null
    mv velvet_* velvet
  else
    echo "Trying to download velvet source code to $RunDir/velvet/version from https://github.com/dzerbino/velvet/trunk"
    cd $RunDir/velvet/version/
    svn co https://github.com/dzerbino/velvet/trunk
    if [ -d $RunDir/velvet/version/trunk ]; then
      mv trunk velvet
    fi
  fi
fi
if [ ! -d $RunDir/velvet/version/velvet ]; then
  echo "Getting velvet source code failed. Please download latest velvet to $RunDir/velvet/ from website: https://www.ebi.ac.uk/~zerbino/velvet/"
  exit 1
fi
cd $RunDir/velvet/version/velvet
make 'CATEGORIES=10' 'MAXKMERLENGTH=301' 'BIGASSEMBLY=1' 'LONGSEQUENCES=1' 'OPENMP=1'
if [ ! -s velveth ] || [ ! -s velvetg ]; then
  echo "Compiling velvet failed"
  exit 1
fi
velvet_version_temp=`./velveth 2>&1 | grep 'Version'`
velvet_version=${velvet_version_temp##* }
if [ "v$velvet_version" = "v" ]; then
  echo "Can not get velvet version. Exit..."
  exit 1
fi
echo "velvet version $velvet_version dettected"
if [ -d  $RunDir/velvet/version/$MachType ]; then 
  rm -rf $RunDir/velvet/version/$MachType
fi
mkdir -p $RunDir/velvet/version/$MachType/bin
cp velvetg velveth $RunDir/velvet/version/$MachType/bin/
cd $RunDir/velvet
rm -rf $RunDir/velvet/version/velvet
mv version v$velvet_version
if [ -d $RunDir/bin ]; then
  mkdir -p $RunDir/bin
fi
ln -sf $RunDir/velvet/v$velvet_version/$MachType/bin/velveth $RunDir/bin/
ln -sf $RunDir/velvet/v$velvet_version/$MachType/bin/velvetg $RunDir/bin/
echo "*******************************************"
echo "************ Important ********************"
echo "*******************************************"
echo "Add to /etc/profile or ~/.bashrc\n\n"
echo "###velvet"
echo "export PATH=$RunDir/velvet/v$velvet_version/$MachType/bin"':$PATH'
echo "\n\n\n"

exit 0
