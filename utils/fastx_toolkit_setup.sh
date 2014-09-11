#!/bin/sh
#Require GNU g++ (>4.1.2), libgtextutils-0.7, PerlIO::gzip, GD::Graph::bars, GNU sed, gnuplot
#GD::Graph::bars requires perl (>5.6.0), GD (>1.19), GD::Text::Align
#GD requires libgd YAML
#APT
#apt-get install gcc g++ pkg-config wget gnuplot
#apt-get install libgd2-xpm-dev build-essential libgd-gd2-perl libconfig-yaml-perl
#YUM
#yum install pkgconfig gcc gcc-c++ wget gnuplot
#yum install gd gd-devel php-gd perl-YAML
RunDir=$(cd `dirname $0`; pwd)
MachType=$(uname -m)

if which fastx_reverse_complement 2>/dev/null; then
  echo "fastx_toolkit already installed"
  exit 0
fi

if [ -d $RunDir/fastx_toolkit ]; then
  mkdir -p $RunDir/fastx_toolkit
fi
cd $RunDir/fastx_toolkit
if [ -d $RunDir/fastx_toolkit/v0.0.14 ]; then
  rm -rf $RunDir/fastx_toolkit/v0.0.14
fi
mkdir -p $RunDir/fastx_toolkit/v0.0.14/src
cd $RunDir/fastx_toolkit/v0.0.14/src
wget https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2
wget https://github.com/agordon/libgtextutils/releases/download/0.7/libgtextutils-0.7.tar.gz
cd $RunDir/fastx_toolkit/v0.0.14/
#install libgtextutils
libgtextutils_exist=0
if [ -s $RunDir/libs/libgtextutils/v0.7/$MachType/lib/libgtextutils-0.7.so.0.0.0 ]; then
  echo "libtextutils-0.7 exists in $RunDir/libs/libgtextutils/v0.7"
  libgtextutils_exist=1
else
  echo "Installing libgtextutils-0.7 ..."
  tar xvf $RunDir/fastx_toolkit/v0.0.14/src/libgtextutils-0.7.tar.gz > /dev/null
  cd libgtextutils-0.7/
  if [ -d $RunDir/libs/libgtextutils/v0.7/$MachType ]; then
    rm -rf $RunDir/libs/libgtextutils/v0.7/$MachType
  fi
  mkdir -p $RunDir/libs/libgtextutils/v0.7/$MachType
  ./configure --prefix=$RunDir/libs/libgtextutils/v0.7/$MachType > configure.log
  make > make.log
  make install > make_install.log
fi
if [ ! -s $RunDir/libs/libgtextutils/v0.7/$MachType/lib/libgtextutils-0.7.so.0.0.0 ]; then
  echo "Compiling libgtextutils-0.7 failed"
  exit 1
fi
rm -rf $RunDir/fastx_toolkit/v0.0.14/libgtextutils-0.7
export CPLUS_INCLUDE_PATH=$RunDir/libs/libgtextutils/v0.7/$MachType/include:$CPLUS_INCLUDE_PATH
export C_INCLUDE_PATH=$RunDir/libs/libgtextutils/v0.7/$MachType/include:$C_INCLUDE_PATH
export LIBRARY_PATH=$RunDir/libs/libgtextutils/v0.7/$MachType/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$RunDir/libs/libgtextutils/v0.7/$MachType/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$RunDir/libs/libgtextutils/v0.7/$MachType/lib/pkgconfig:$PKG_CONFIG_PATH
###install fastxtoolkit
echo "Installing fastx_toolkit ..."
cd $RunDir/fastx_toolkit/v0.0.14/
if [ -d $RunDir/fastx_toolkit/v0.0.14/fastx_toolkit-0.0.14 ]; then
  rm -rf $RunDir/fastx_toolkit/v0.0.14/fastx_toolkit-0.0.14
fi
tar xvf $RunDir/fastx_toolkit/v0.0.14/src/fastx_toolkit-0.0.14.tar.bz2 > /dev/null
cd $RunDir/fastx_toolkit/v0.0.14/fastx_toolkit-0.0.14
./configure --prefix=$RunDir/fastx_toolkit/v0.0.14/$MachType > configure.log
make > make.log
make install > make_install.log
if [ ! -s $RunDir/fastx_toolkit/v0.0.14/$MachType/bin/fastx_reverse_complement ]; then
  echo "Compiling fastx_toolkit failed"
  exit 1
fi
cd $RunDir/fastx_toolkit/v0.0.14/
rm -rf $RunDir/fastx_toolkit/v0.0.14/fastx_toolkit-0.0.14
if [ -d $RunDir/bin ]; then
  mkdir -p $RunDir/bin
fi
ln -sf $RunDir/fastx_toolkit/v0.0.14/$MachType/bin/fastx_reverse_complement $RunDir/bin/
ln -sf $RunDir/fastx_toolkit/v0.0.14/$MachType/bin/fastx_trimmer $RunDir/bin/
echo "*******************************************"
echo "************ Important ********************"
echo "*******************************************"
echo "Add to /etc/profile or ~/.bashrc\n\n"
if [ $libgtextutils_exist = 0 ]; then
  echo "###libgtextutils-0.7"
  echo "export CPLUS_INCLUDE_PATH=$RunDir/libs/libgtextutils/v0.7/$MachType/include"':$CPLUS_INCLUDE_PATH'
  echo "export C_INCLUDE_PATH=$RunDir/libs/libgtextutils/v0.7/$MachType/include"':$C_INCLUDE_PATH'
  echo "export LIBRARY_PATH=$RunDir/libs/libgtextutils/v0.7/$MachType/lib"':$LIBRARY_PATH'
  echo "export LD_LIBRARY_PATH=$RunDir/libs/libgtextutils/v0.7/$MachType/lib"':$LD_LIBRARY_PATH'
  echo "export PKG_CONFIG_PATH=$RunDir/libs/libgtextutils/v0.7/$MachType/lib/pkgconfig"':$PKG_CONFIG_PATH'
fi
echo "###fastx_toolkit"
echo "export PATH=$RunDir/fastx_toolkit/v0.0.14/$MachType/bin"':$PATH'
echo "\n\n\n"

exit 0
