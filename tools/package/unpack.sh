#!/bin/bash

# usage: ./unpack.sh /path/to/public/repo
SOURCE=$1
if [ -z $SOURCE ]; then
  echo "must input source"
  exit
fi

# clone the repository
CLONE="/tmp/feasst"
if [ -d $CLONE ]; then
  echo "clone directory must not previously exist"
  exit
fi
if [ -z $CLONE ]; then
  echo "CLONE directory must not be empty "$CLONE
  exit
fi
git clone $SOURCE $CLONE

# strip git from the clone
rm -rf $CLONE/.git

# strip stubs from the clone
filelist=""
if [ -d stubs ]; then
  pushd stubs
  filelist=`find . -type l`
  popd
fi
for file in $filelist; do
  rm $CLONE/$file
done

# remove private files from the clone
while IFS= read -r var
do
  rm -r $CLONE/${var}*
done < "private_files.txt"

# copy the repo to destination
shopt -s dotglob nullglob   #this command moves hidden files
rsync -r $CLONE/* ../../
#mv $CLONE/* ../../

echo "##################################################"
echo "summarize differences between clone and desination"
echo "##################################################"

diff -qNr $CLONE $SOURCE | grep -v "/.git/"

# clean the clone
rm -rf $CLONE

