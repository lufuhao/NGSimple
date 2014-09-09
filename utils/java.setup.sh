#!/bin/sh
###configure
###download JDK according x86 or x64 to $RunDir/java/
#http://www.oracle.com/technetwork/java/javase/downloads/index.html
###configure end

if which java 2>/dev/null; then
  echo "JAVA already installed"
  exit 1
fi

RunDir=$(cd `dirname $0`; pwd)
MachType=$(uname -m)
if [ -d $RunDir/java ]; then
  echo "Dir exists: $RunDir/java"
else
  mkdir -p $RunDir/java
fi
cd $RunDir/java
if [ -s jdk*-linux-*.tar.gz ]; then
  jdkpackage=$(ls jdk*-linux-*.tar.gz)
  echo "JDK file detected: $jdkpackage"
else
  echo "Please download JDK package to $RunDir/java/ from website http://www.oracle.com"
  exit 0
fi

jdk_basename=${jdkpackage%%.*}
jdk_MachType=${jdk_basename##*-}
echo "JDK_machtype: $jdk_MachType"
jdk_version1=${jdk_basename%-linux*}
jdk_version2=${jdk_version1#*-}
jdk_major=${jdk_version2%u*}
jdk_minor=${jdk_version2#*u}
jdk_version="1.$jdk_major.0.$jdk_minor"
echo "jdk version: $jdk_version"
mkdir -p $RunDir/java/v$jdk_version/src
mv $jdkpackage $RunDir/java/v$jdk_version/src/
cd $RunDir/java/v$jdk_version
tar zxvf src/$jdkpackage >/dev/null
mv jdk* $jdk_MachType
if [ -d $RunDir/bin ]; then
  echo "java is soft-linked to $RunDir/bin"
else
  mkdir -p $RunDir/bin
fi
cd $RunDir/java/v$jdk_version/$jdk_MachType/bin
ln -sf $RunDir/java/v$jdk_version/$jdk_MachType/bin/java $RunDir/bin/
echo "*******************************************"
echo "************ Important ********************"
echo "*******************************************"
echo "Add to /etc/profile or ~/.bashrc\n\n"
echo "###JAVA"
echo "export JAVA_HOME=$RunDir/java/v$jdk_version/$jdk_MachType"
echo 'export JRE_HOME=$JAVA_HOME/jre'
echo 'export PATH=$JAVA_HOME/bin:$JRE_HOME/bin:$PATH'
echo 'export CLASSPATH=.:$JAVA_HOME/lib:$JRE_HOME/lib:$CLASSPATH'
echo "\n\n\n"
exit 1
