#!/bin/sh
#Require g++, zlib, GNU ncurses library, subversion (for bamtools), cmake (>=2.6.4), sparsehash, gsl
#Sparsehash require: aclocal-1.11, automake-1.11, autoconf, autoheader
#APT
#apt-get install zlib1g zlib1g-dev
#apt-get install libncurses5 libncurses5-dev
#apt-get install sparsehash libsparsehash-dev
#apt-get install libgsl0ldbl libgsl0-dev
#apt-get install automake
#YUM
#yum install zlib zlib-devel
#yum install ncurses ncurses-devel
#yum install sparsehash sparsehash-devel
#yum install gsl gsl-devel
#?
RunDir=$(cd `dirname $0`; pwd)
MachType=$(uname -m)

if which fastq-join 2>/dev/null; then
  echo "ea-utils already installed"
  exit 0
fi

#bamtools
bamtools_exist=0
if which bamtools 2>/dev/null; then
  echo "Bamtools exists, which is required to build sam-stats in ea-utils package"
  bamtools_exist=1
else
  echo "Installing required bamtools"
  if [ ! -d $RunDir/bamtools ]; then
    mkdir -p $RunDir/bamtools
  fi
  cd $RunDir/bamtools
  if [ -d $RunDir/bamtools/bamtools ]; then
    rm -rf $RunDir/bamtools/bamtools
  fi
  mkdir -p $RunDir/bamtools/bamtools
  cd $RunDir/bamtools/bamtools
  svn co https://github.com/pezmaster31/bamtools/trunk > /dev/null
  if [ -d $RunDir/bamtools/bamtools/trunk ]; then
    mv trunk $MachType
  else
    echo "Can not find BamTools source trunk in $RunDir/bamtools/bamtools/branches"
    exit 1
  fi
  cd $RunDir/bamtools/bamtools/$MachType
  if [ -d $RunDir/bamtools/bamtools/$MachType/build ]; then
    rm -rf $RunDir/bamtools/bamtools/$MachType/build
  fi
  mkdir -p $RunDir/bamtools/bamtools/$MachType/build
  cd $RunDir/bamtools/bamtools/$MachType/build
  cmake .. > $RunDir/bamtools/bamtools/$MachType/cmake.logs
  make > $RunDir/bamtools/bamtools/$MachType/make.logs
  cd $RunDir/bamtools/bamtools/$MachType
  if [ ! -s $RunDir/bamtools/bamtools/$MachType/bin/bamtools-* ]; then
    echo "Compiling BamTools failed"
    exit 1
  fi
  cd $RunDir/bamtools/bamtools/$MachType/bin
  bamtools_bin=`ls bamtools-*`
  bamtools_version=${bamtools_bin##*-}
  if [ "v$bamtools_version" = "v" ]; then
    echo "Can not guess bamtools version"
    exit 1
  fi
  cd $RunDir/bamtools/
  if [ -d $RunDir/bamtools/v$bamtools_version ]; then
    rm -rf $RunDir/bamtools/v$bamtools_version
  fi
  mv bamtools v$bamtools_version
  if [ -d $RunDir/bin ]; then
    mkdir -p $RunDir/bin
  fi
  ln -sf $RunDir/bamtools/v$bamtools_version/$MachType/bin/$bamtools_bin $RunDir/bin/bamtools
  export PATH=$RunDir/bamtools/v$bamtools_version/$MachType/bin:$PATH
  export CPLUS_INCLUDE_PATH=$RunDir/bamtools/v$bamtools_version/$MachType/include:$CPLUS_INCLUDE_PATH
  export C_INCLUDE_PATH=$RunDir/bamtools/v$bamtools_version/$MachType/include:$C_INCLUDE_PATH
  export LIBRARY_PATH=$RunDir/bamtools/v$bamtools_version/$MachType/lib:$LIBRARY_PATH
  export LD_LIBRARY_PATH=$RunDir/bamtools/v$bamtools_version/$MachType/lib:$LD_LIBRARY_PATH
fi
###
###
###ea-utils
sparsehash_exist=0
if [ ! -d $RunDir/eautils ]; then
  mkdir -p $RunDir/eautils
fi
cd $RunDir/eautils
if [ ! -s ea-utils*.tar.gz ]; then
  echo "Please download the EAutils source code to $RunDir/eautils/ from website: https://code.google.com/p/ea-utils/"
  exit 1
fi
eautils_package=`ls ea-utils*.tar.gz`
eautils_basename=${eautils_package%.*.*}
eautils_version=${eautils_basename#*.}
if [ -d $RunDir/eautils/v$eautils_version ]; then
  rm -rf $RunDir/eautils/v$eautils_version
fi
mkdir -p $RunDir/eautils/v$eautils_version/src
mv $eautils_package $RunDir/eautils/v$eautils_version/src/
cd $RunDir/eautils/v$eautils_version/
tar xvf $RunDir/eautils/v$eautils_version/src/$eautils_package > /dev/null
mv ea-utils* sourcecode
cd $RunDir/eautils/v$eautils_version/sourcecode/
if [ -s $RunDir/libs/sparsehash/v2.0.2/x86_64/include/sparsehash/sparse_hash_map ]; then
  sparsehash_exist=1
  echo "Google sparseharsh v2.0.2 exists in $RunDir/libs/sparsehash/v2.0.2/x86_64"
else
  echo "Installing Google Sparsehash..."
  mkdir -p $RunDir/libs/sparsehash/v2.0.2/x86_64
  cd $RunDir/eautils/v$eautils_version/sourcecode/sparsehash-2.0.2
  ./configure --prefix=$RunDir/libs/sparsehash/v2.0.2/x86_64 > configure.log
  make > make.log
  make install > make_install.log
  if [ ! -s $RunDir/libs/sparsehash/v2.0.2/x86_64/include/sparsehash/sparse_hash_map ]; then
    echo "Compiling Sparsehash failed"
    exit 1
  fi

fi
export CPLUS_INCLUDE_PATH=$RunDir/libs/sparsehash/v2.0.2/x86_64/include:$CPLUS_INCLUDE_PATH
export C_INCLUDE_PATH=$RunDir/libs/sparsehash/v2.0.2/x86_64/include:$C_INCLUDE_PATH
export PKG_CONFIG_PATH=$RunDir/libs/sparsehash/v2.0.2/x86_64/lib/pkgconfig:$PKGCONFIG
cd $RunDir/eautils/v$eautils_version/sourcecode/
echo "Installing EA-utils..."
PREFIX=$RunDir/eautils/v$eautils_version/$MachType make install > make.log
if [ ! -s $RunDir/eautils/v$eautils_version/$MachType/bin/fastq-join ]; then
  echo "Compiling ea-utils failed"
  exit 1
fi
cd $RunDir/eautils/v$eautils_version/
rm -rf $RunDir/eautils/v$eautils_version/sourcecode
if [ -d $RunDir/bin ]; then
  mkdir -p $RunDir/bin
fi
ln -sf $RunDir/eautils/v$eautils_version/$MachType/bin/fastq-join $RunDir/bin/
echo "*******************************************"
echo "************ Important ********************"
echo "*******************************************"
echo "Add to /etc/profile or ~/.bashrc\n\n"
if [ $sparsehash_exist = 0 ]; then
  echo "###Sparsehash"
  echo "export CPLUS_INCLUDE_PATH=$RunDir/libs/sparsehash/v2.0.2/x86_64/include"':$CPLUS_INCLUDE_PATH'
  echo "export C_INCLUDE_PATH=$RunDir/libs/sparsehash/v2.0.2/x86_64/include"':$C_INCLUDE_PATH'
  echo "export PKG_CONFIG_PATH=$RunDir/libs/sparsehash/v2.0.2/x86_64/lib/pkgconfig"':$PKGCONFIG'
fi
if [ $bamtools_exist = 0 ]; then
  echo "###BamTools"
  echo "export PATH=$RunDir/bamtools/v$bamtools_version/$MachType/bin"':$PATH'
  echo "export CPLUS_INCLUDE_PATH=$RunDir/bamtools/v$bamtools_version/$MachType/include"':$CPLUS_INCLUDE_PATH'
  echo "export C_INCLUDE_PATH=$RunDir/bamtools/v$bamtools_version/$MachType/include"':$C_INCLUDE_PATH'
  echo "export LIBRARY_PATH=$RunDir/bamtools/v$bamtools_version/$MachType/lib"':$LIBRARY_PATH'
  echo "export LD_LIBRARY_PATH=$RunDir/bamtools/v$bamtools_version/$MachType/lib"':$LD_LIBRARY_PATH'
fi
echo "###ea-utils"
echo "export PATH=$RunDir/eautils/v$eautils_version/$MachType/bin"':$PATH'
echo "\n\n\n"

exit 0
