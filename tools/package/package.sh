#!/bin/bash

# usage: ./package.sh /path/to/public/repo
DEST=$1
if [ -z $DEST ]; then
  echo "must input destination"
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
git clone ../../ $CLONE

# strip git from the clone
rm -rf $CLONE/.git

# copy stubs to the clone
# -L flag copies actual files from soft links
cp -rL stubs/* $CLONE/

# remove private files from the clone
while IFS= read -r var
do
  rm -r $CLONE/${var}*
done < "private_files.txt"

# # copy the repo to destination
# shopt -s dotglob nullglob   #this command moves hidden files
# mv $CLONE $DEST

# copy the repo to the destination
rsync -azv $CLONE/ $DEST/

echo "##################################################"
echo "summarize differences between clone and desination"
echo "##################################################"

diff -qNr $CLONE $DEST | grep -v "/.git/"

# clean the clone in preparation of next package
rm -rf $CLONE

