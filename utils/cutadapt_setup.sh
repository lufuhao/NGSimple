#!/bin/sh
#Require python (2.6<=version<3), c compiler, python-dev, Cython
###APT
#apt-get install python-dev cython
###YUM
#yum install python-devel cython
RunDir=$(cd `dirname $0`; pwd)
MachType=$(uname -m)

if which cutadapt 2>/dev/null; then
  echo "Cutadapt already installed"
  exit 1
fi

if [ ! -d $RunDir/cutadapt ]; then
  mkdir -p $RunDir/cutadapt
fi
cd $RunDir/cutadapt
if [ ! -s cutadapt-*.tar.gz ]; then
  wget https://pypi.python.org/packages/source/c/cutadapt/cutadapt-1.4.1.tar.gz
fi
if [ -s cutadapt-*.tar.gz ]; then
  cutadapt_package=`ls cutadapt-*.tar.gz`
  cutadapt_basename=${cutadapt_package%.*.*}
  cutadapt_version=${cutadapt_basename##*-}
  if [ "v$cutadapt_version" = "v" ]; then
    echo "Can not get cutadapt version, try to continue..."
  else
    echo "cutadapt version $cutadapt_version detected"
  fi
else
  echo "Downloading CutAdapt failed"
  exit 0
fi
if [ -d $RunDir/cutadapt/$cutadapt_basename ]; then
 rm -rf $RunDir/cutadapt/$cutadapt_basename
fi
tar xvf $cutadapt_package > /dev/null
cd $RunDir/cutadapt/$cutadapt_basename
python setup.py build > build.log
log1=`python setup.py install --user | grep "changing mode of "`
log3=${log1#changing mode of }
log4=${log3%'/cutadapt to 775'}
echo "cutadapt was installed to: $log4"
if [ -d $RunDir/bin ]; then
  mkdir -p $RunDir/bin
fi
if [ -s $log4/cutadapt ]; then
  echo "Detecting cutadapt path ...1..."
  ln -sf $log4/cutadapt $RunDir/bin/
  cutadapt_path=$log4
elif [ -s $HOME/.local/bin/cutadapt ]; then
  echo "Detecting cutadapt path ...2..."
  ln -sf $HOME/.local/bin/cutadapt $RunDir/bin/
  cutadapt_path="$HOME/.local/bin"
else
  echo "Can not find installed cutadapt path"
  exit 0
fi
echo "*******************************************"
echo "************ Important ********************"
echo "*******************************************"
echo "Add to /etc/profile or ~/.bashrc\n\n"
echo "###cutadapt"
echo "export PATH=$cutadapt_path"':$PATH'
echo "\n\n\n"
