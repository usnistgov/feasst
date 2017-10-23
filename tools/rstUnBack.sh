#!/usr/bin/env bash
files=`ls * | grep -v bak`
for file in $files; do
  #echo "$file $file.bak"
  mv ${file}.bak $file
done
