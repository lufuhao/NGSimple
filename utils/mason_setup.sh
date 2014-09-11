#!/bin/sh
#1. Mason
#Require: subversion, cmake, make, C++ compiler, java, zlib, bzip2, boost
#Optional: LEMON library,
#APT
#apt-get install subversion
#apt-get install g++
#apt-get install cmake
###zlib
#apt-get install zlib1g zlib1g-dev
###boost
#apt-get install libboost-dev
###bzip2 lib
#apt-get install libbz2-dev
###LEMON library
#apt-get install lemon

#YUM
#yum install gcc gcc-c++
#yum install subversion
#yum install cmake
###zlib
#yum install zlib zlib-devel
###boost
#yum install boost boost-devel boost-doc
###bzip2 lib
#yum install bzip2 bzip2-devel

if which mason 2>/dev/null; then
  echo "Mason already installed"
  exit 0
fi

RunDir=$(cd `dirname $0`; pwd)
MachType=$(uname -m)
if [ -d $RunDir/mason ]; then
  echo "Dir exists: $RunDir/mason"
else
  mkdir -p $RunDir/mason
fi

cd $RunDir/mason
if [ -d $RunDir/mason/seqan-trunk ]; then
  echo "Mason source code detected in folder $RunDir/mason/seqan-trunk"
else
  if [ -s mason*.tar.gz ]; then
    MasonTarGz=`ls mason*.tar.gz`
    echo "Un-compress file: $MasonTarGz"
    tar zxvf $MasonTarGz > /dev/null
  else
    echo "Downloading Mason using subversion"
    svn co http://svn.seqan.de/seqan/trunk seqan-trunk > /dev/null
  fi
  if [ -d $RunDir/mason/seqan-trunk ]; then
    echo "Mason source code detected in folder $RunDir/mason/seqan-trunk"
  else
    echo "Can not find Mason source code"
    exit 1
  fi
fi
if [ -s mason_v2010.tar.gz ]; then
  echo "mason backup file exists: mason_v2010.tar.gz"
else
  tar zcvf mason_v2010.tar.gz seqan-trunk/ > /dev/null
  echo "Mason source code was compressed to file: mason_v2010.tar.gz"
fi
if [ -d $RunDir/mason/v2010 ]; then
  rm -rf $RunDir/mason/v2010
fi
mkdir -p $RunDir/mason/v2010/$MachType
cd $RunDir/mason/v2010/$MachType
rm -rf *
cmake $RunDir/mason/seqan-trunk > cmake.log
make mason > make.log
if [ -s $RunDir/mason/v2010/$MachType/bin/mason ]; then
  ln -sf $RunDir/mason/v2010/$MachType/bin/mason $RunDir/bin/
else
  echo "Compiling Mason failed"
  exit 1
fi
mv $RunDir/mason/seqan-trunk/core/include $RunDir/mason/v2010/$MachType/
rm -rf $RunDir/mason/seqan-trunk
echo "*******************************************"
echo "************ Important ********************"
echo "*******************************************"
echo "Add to /etc/profile or ~/.bashrc\n\n"
echo "###Mason"
echo "export PATH=$RunDir/mason/v2010/$MachType/bin"':$PATH'
echo "export CPLUS_INCLUDE_PATH=$RunDir/mason/v2010/$MachType/include"':$CPLUS_INCLUDE_PATH'
echo "\n\n\n"
exit 0
