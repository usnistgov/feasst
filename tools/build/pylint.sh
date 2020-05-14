file=$(basename $0 .sh)
echo "*** SRC ***" > $file.log 2>&1
pylint ../src/*.py >> $file.log 2>&1
files=`find ../tutorial/ -name "test.py"`
for pyfile in $files; do
  echo "*** $pyfile ***" >> $file.log 2>&1
  pylint $pyfile >> $file.log 2>&1
done
# find all py files in tutorial
echo "**** $file ****" >> summary.log
tail -1 $file.log >> summary.log
echo "" >> summary.log
