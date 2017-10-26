#!/usr/bin/env bash
# This script makes the feasst.i and factories.cc files, depending on which header files are present.
# usage: ./makeFactory.sh

# Note that the following classes must be defined in a very specific order
factoryBaseClasses="pair criteria analyze trial random"
baseClasses="custom_exception random base table histogram accumulator space pair barrier shape criteria analyze trial mc"
nonClasses="functions ui_abbreviated"

##############################
echo "Generating feasst.i"
##############################

swigfile="feasst.i"
cat << EOF > $swigfile
/* NOTE: This file is created by tools/makeFactory.sh. Modify this file instead.
   See buildtemplate/feasst.i.example for an example file. */

%module feasst

%ignore CriteriaWLTMMC::lnPIrwsatwrap;
%ignore CriteriaWLTMMC::lnPIrwnmxwrap;
%ignore MC::boyleminwrap;

%{
EOF

classes="$baseClasses $nonClasses"
printfile=$swigfile

prepend="#include \""
append="\""
printHeaders()
{
  for class in $classes; do
    for header in `ls -vd ../src/${class}*.h 2> /dev/null`; do
      echo ${prepend}`basename ${header}`${append} >> $printfile
    done
  done
}

printHeaders

cat << EOF >> $swigfile
#ifdef XDRFILE_H_
  extern "C" {
    #include "xdrfile.h"
    #include "xdrfile_xtc.h"
    #include "xdrfile_trr.h"
  }
#endif // XDRFILE_H_
%}

%include "std_vector.i"
%include "std_shared_ptr.i"
%include "std_string.i"
%include "std_iostream.i"
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
using namespace std;
// using std::vector;

%pythonnondynamic;

%shared_ptr(std::ifstream);
namespace std{
  class ifstream {
  };
}

%shared_ptr(Base);
EOF

printTemplates()
{
  for class in $classes; do
    for header in `ls -vd ../src/${class}*.h 2> /dev/null`; do
      classNameList=`grep -h "class \(.*\) : public \(.*\) {" $header | sed 's/class \(.*\) : public \(.*\) {/\1/'`
      for className in $classNameList; do
        echo "%shared_ptr("$className");" >> $printfile
      done
    done
  done
}

prepend=""
append=""
printTemplates

prepend="%include "
append=""
printHeaders

##############################
echo "Generating factories.cc"
##############################

ccfile="../src/factories.cc"
cat << EOF > $ccfile
/* Factory method for Pair, Criteria, Trial and Analyze.
   This file is created by tools/makeFactory.sh. Modify this file instead.
   See src/factories.cc.example for an example file. */

EOF
printfile=$ccfile
prepend="#include \""
append="\""
for classes in $factoryBaseClasses; do
  printHeaders
  echo "" >> $ccfile
done

cat << EOF >> $ccfile
#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Pair* makePair(Space* space, const char* fileName) {
  ASSERT(fileExists(fileName),
    "restart file(" << fileName << ") doesn't exist");
  string typestr = fstos("className", fileName);
  Pair* pair = NULL;
  if (1 == 0) { // filler for a loop of "} else if {"
EOF

printFactory()
{
  for header in `ls -vd ../src/${class}_*.h 2> /dev/null`; do
    className=`grep "class " $header | grep "{" | head -1 | awk '{print $2}'`
    if [ $shrptr == "no" ]; then
      echo "  } else if (typestr.compare(\""$className"\") == 0) { "$class" = new "${className}${args}";" >> $printfile
    else
      echo "  } else if (typestr.compare(\""$className"\") == 0) { "$class" = make_shared<"${className}">"${args}";" >> $printfile
    fi
  done
}
appendPrint()
{
  echo "  } else { ASSERT(0, \"unrecognized "$class"(\" << typestr << \") in factory\"); }" >> $printfile
  echo "  return "$class";" >> $printfile
  echo "}" >> $printfile
  echo "" >> $printfile
}

class=pair
args="(space, fileName)"
shrptr=no
printFactory
appendPrint

cat << EOF >> $ccfile
Criteria* makeCriteria(const char* fileName) {
  ASSERT(fileExists(fileName),
    "restart file(" << fileName << ") doesn't exist");
  string typestr = fstos("className", fileName);
  Criteria* criteria = NULL;
  if (typestr.compare("Criteria") == 0) {
    criteria = new CriteriaMetropolis(fileName);
EOF

class=criteria
args="(fileName)"
shrptr=no
printFactory
appendPrint

cat << EOF >> $ccfile
shared_ptr<Trial> makeTrial(Pair* pair, Criteria* criteria,
  const char* fileName) {
  ASSERT(fileExists(fileName),
         "restart file(" << fileName << ") doesn't exist");
  string typestr = fstos("className", fileName);
  shared_ptr<Trial> trial;
  if (1 == 0) { // filler for a loop of "} else if {"
EOF

class=trial
shrptr=yes
args="(fileName, pair, criteria)"
printFactory
appendPrint

cat << EOF >> $ccfile
shared_ptr<Analyze> makeAnalyze(Pair* pair,
  const char* fileName) {
  ASSERT(fileExists(fileName),
         "restart file(" << fileName << ") doesn't exist");
  string typestr = fstos("className", fileName);
  shared_ptr<Analyze> analyze;
  if (typestr.compare("Analyze") == 0) {
    analyze = make_shared<Analyze>(pair, fileName);
EOF

class=analyze
args="(pair, fileName)"
shrptr=yes
printFactory
appendPrint

cat << EOF >> $ccfile
#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_
EOF

#####################################
echo "Generating documentation *.rst"
#####################################

printRSTHeader()
{
  className=`grep "class " $header | grep "{" | head -1 | awk '{print $2}'`
  printfile=${className}.rst
  echo "   $className" >> $tocfile
  echo "$className" > $printfile
  echo "=====================================================" >> $printfile
  echo "" >> $printfile
  if [ -z "$nonClassFlag" ]; then
    echo ".. doxygenclass:: $className" >> $printfile
    echo "   :project: FEASST" >> $printfile
    echo "   :members:" >> $printfile
  else
    echo ".. doxygenfile:: `basename $header`" >> $printfile
    echo "   :project: FEASST" >> $printfile
  fi
}

printDoc()
{
  arr=(`ls -vd ../src/${class}_*.h 2> /dev/null`)
  # if array finds derived classes, set up toc
  header="../src/${class}.h"
  className=`grep "class " $header | grep "{" | head -1 | awk '{print $2}'`
  tocfile=${class}toc.rst
  echo "$className Classes" > $tocfile
  echo "=======================================================" >> $tocfile
  echo "" >> $tocfile
  echo ".. toctree::" >> $tocfile
  printRSTHeader
  if (( ${#arr[@]} > 0 )); then
    for header in ${arr[@]}; do
      printRSTHeader
    done
    echo "   ${class}toc" >> $apifile
  else
    echo "   $className" >> $apifile
  fi
}

# begin printing the api file, which needs to know whether to include one
# class, or a toc for may derived classes
apifile=api.rst
cat <<-EOF > $apifile
***********************************
Application program interface (API)
***********************************

This section describes the API which is available for use via C++ or python.

.. toctree::
EOF

# print the rst for each class and also add to the api
nonClassFlag=""
for class in $baseClasses; do
  printDoc
done
nonClassFlag="True"
for class in $nonClasses; do
  printDoc
done

