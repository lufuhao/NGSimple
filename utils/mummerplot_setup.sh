#!/bin/sh
#Require: make (GNU make 3.79.1), perl (PERL 5.6.0), sh (GNU sh 1.14.7), csh (tcsh 6.10.00), g++ (GNU gcc 2.95.3), sed (GNU sed 3.02), awk (GNU awk 3.0.4), ar (GNU ar 2.9.5)
#Optional: fig2dev (fig2dev 3.2.3), gnuplot (gnuplot 4.0), xfig (xfig 3.2)
##APT
#apt-get install tcsh xfig gnuplot
###YUM
#yum install tcsh xfig gnuplot
RunDir=$(cd `dirname $0`; pwd)
MachType=$(uname -m)

if which mummerplot 2>/dev/null; then
  echo "Mummerplot already installed"
  exit 0
fi

if [ ! -d $RunDir/mummer ]; then
  mkdir -p $RunDir/mummer
fi
cd $RunDir/mummer
if [ ! -s MUMmer*.tar.gz ]; then
  echo "Trying to download mummer to $RunDir/mummer/ from website: http://sourceforge.net/projects/mummer"
  wget http://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz
fi
if [ ! -s MUMmer*.tar.gz ]; then
  echo "Mummer source code is not present and can not bed downloaded"
  exit 1
fi
mummer_package=`ls MUMmer*.tar.gz`
mummer_basename=${mummer_package%.*.*}
mummer_version=${mummer_basename#MUMmer}
if [ "v$mummer_version" = "v" ]; then
  echo "Can not detect MUMmer version"
  exit 1
fi
echo "Mummer version $mummer_version detected"
if [ -d $RunDir/mummer/v$mummer_version ]; then
  rm -rf $RunDir/mummer/v$mummer_version
fi
mkdir -p $RunDir/mummer/v$mummer_version/src
mv $mummer_package $RunDir/mummer/v$mummer_version/src/
cd $RunDir/mummer/v$mummer_version/
tar xvf $RunDir/mummer/v$mummer_version/src/$mummer_package >/dev/null
if [ ! -d MUMmer* ]; then
  echo "Failed to un-compress $RunDir/mummer/v$mummer_version/src/$mummer_package"
  exit 1
fi
mv MUMmer* $MachType
cd $RunDir/mummer/v$mummer_version/$MachType
echo "Installing MUMmer ..."
make check >make.log
make install >makeinstall.log
if [ ! -s nucmer ] || [ ! -s delta-filter ] || [ ! -s mummerplot ]; then
  echo "Compiling MUMmer failed"
  exit 1
fi
if [ -d $RunDir/bin ]; then
  mkdir -p $RunDir/bin
fi
ln -sf $RunDir/mummer/v$mummer_version/$MachType/nucmer $RunDir/bin/
ln -sf $RunDir/mummer/v$mummer_version/$MachType/delta-filter $RunDir/bin/
ln -sf $RunDir/mummer/v$mummer_version/$MachType/mummerplot $RunDir/bin/
echo "*******************************************"
echo "************ Important ********************"
echo "*******************************************"
echo "Add to /etc/profile or ~/.bashrc\n\n"
echo "###Mummer"
echo "export PATH=$RunDir/mummer/v$mummer_version/$MachType"':$PATH'
echo "\n\n\n"
exit 0
