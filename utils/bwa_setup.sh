#!/bin/sh
#Require zlib
RunDir=$(cd `dirname $(readlink -f $0)`; pwd)
MachType=$(uname -m)

if which bwa >/dev/null; then
  echo "BWA fund installed"`which bwa`
else
  echo "BWA NOT found"
  exit 1
fi

if [ ! -d $RunDir/bwa ]; then
mkdir -p $RunDir/bwa
fi
cd $RunDir/bwa
if [ -s bwa-*.tar.bz2 ]; then
  BWA_package=`ls bwa-*.tar.bz2`
  echo "BWA source code detected: $BWA_package"
  BWA_basename=${BWA_package%.*.*}
  BWA_version=${BWA_basename#*-}
  if [ -d $RunDir/bwa/v$BWA_version ]; then
    rm -rf $RunDir/bwa/v$BWA_version
  fi
  mkdir -p $RunDir/bwa/v$BWA_version/src
  mv $RunDir/bwa/$BWA_package $RunDir/bwa/v$BWA_version/src/
  cd $RunDir/bwa/v$BWA_version
  tar xvf $RunDir/bwa/v$BWA_version/src/$BWA_package > /dev/null
  mv bwa* $MachType
  cd $RunDir/bwa/v$BWA_version/$MachType
  make > make.log
else
  echo "Trying to download BWA from GitHub"
  if [ -d $RunDir/bwa/trunk ]; then
    rm -rf $RunDir/bwa/trunk
  fi
  svn co https://github.com/lh3/bwa/trunk > /dev/null
  if [ -d $RunDir/bwa/trunk ]; then
    if [ -d $RunDir/bwa/version ]; then
      rm -rf $RunDir/bwa/version
    fi
    mkdir version
    mv trunk version/$MachType
    cd $RunDir/bwa/version/$MachType
    make > make.log
    if [ ! -s $RunDir/bwa/version/$MachType/bwa ]; then
      echo "Compiling BWA failed. Exit..."
      exit 1
    fi
    BWA_version_temp=`./bwa 2>&1 | grep 'Version'`
    BWA_version=${BWA_version_temp##* }
    if [ "v$BWA_version" = "v" ]; then
      echo "Can not get BWA version. Exit..."
      exit 1
    fi
    cd $RunDir/bwa
    mv version v$BWA_version
    cd $RunDir/bwa/v$BWA_version/$MachType
  else
    echo "Can not download BWA using subversion. Please manually download BWA source code to $RunDir/bwa from http://sourceforge.net/projects/bio-bwa"
    exit 1
  fi
fi
echo "BWA version: $BWA_version"
mkdir -p 'bin' 'include' 'lib' 'man/man1'
cp $RunDir/bwa/v$BWA_version/$MachType/bwa ./bin/
cp $RunDir/bwa/v$BWA_version/$MachType/qualfa2fq.pl ./bin/
cp $RunDir/bwa/v$BWA_version/$MachType/xa2multi.pl ./bin/
cp *.h ./include/
cp $RunDir/bwa/v$BWA_version/$MachType/libbwa.a ./lib/
cp $RunDir/bwa/v$BWA_version/$MachType/bwa.1 ./man/man1/
ln -sf $RunDir/bwa/v$BWA_version/$MachType/bin/bwa $RunDir/bin/
echo "*******************************************"
echo "************ Important ********************"
echo "*******************************************"
echo "Add to /etc/profile or ~/.bashrc\n\n"
echo "###BWA"
echo "export PATH=$RunDir/bwa/v$BWA_version/$MachType/bin"':$PATH'
echo "export CPLUS_INCLUDE_PATH=$RunDir/bwa/v$BWA_version/$MachType/include"':$CPLUS_INCLUDE_PATH'
echo "export C_INCLUDE_PATH=$RunDir/bwa/v$BWA_version/$MachType/include"':$C_INCLUDE_PATH'
echo "export LIBRARY_PATH=$RunDir/bwa/v$BWA_version/$MachType/lib"':$LIBRARY_PATH'
echo "export MAN_PATH=$RunDir/bwa/v$BWA_version/$MachType/man"':$MAN_PATH'
echo "\n\n\n"

exit 0
