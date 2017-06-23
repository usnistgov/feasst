#!/bin/bash
# This script makes the feasst.i and factories.cc files, depending on which header files are present.
# usage: ./makeFactory.sh

# Note that the following classes must be defined in a very specific order
factoryBaseClasses="pair criteria analyze trial"
baseClasses="functions custom_exception base table histogram accumulator space pair barrier shape criteria analyze trial mc mc_wltmmc ui_abbreviated"

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

classes=$baseClasses
printfile=$swigfile

prepend="#include \""
append="\""
printHeaders()
{
  for class in $classes; do
    for header in `ls -vd ../src/${class}*.h`; do
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

EOF

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

Pair* Pair::makePair(Space* space, const char* fileName) {
  ASSERT(fileExists(fileName),
    "restart file(" << fileName << ") doesn't exist");
  string typestr = fstos("className", fileName);
  Pair* pair = NULL;
  if (1 == 0) { // filler for a loop of "} else if {"
EOF

printFactory()
{
  for header in `ls -vd ../src/${class}_*.h`; do
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
Criteria* Criteria::makeCriteria(const char* fileName) {
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
shared_ptr<Trial> Trial::makeTrial(Space* space, Pair* pair, Criteria* criteria,
  const char* fileName) {
  ASSERT(fileExists(fileName),
         "restart file(" << fileName << ") doesn't exist");
  string typestr = fstos("className", fileName);
  shared_ptr<Trial> trial;
  if (1 == 0) { // filler for a loop of "} else if {"
EOF

class=trial
shrptr=yes
args="(fileName, space, pair, criteria)"
printFactory
appendPrint

cat << EOF >> $ccfile
shared_ptr<Analyze> Analyze::makeAnalyze(Space* space, Pair* pair,
  const char* fileName) {
  ASSERT(fileExists(fileName),
         "restart file(" << fileName << ") doesn't exist");
  string typestr = fstos("className", fileName);
  shared_ptr<Analyze> analyze;
  if (typestr.compare("Analyze") == 0) {
    analyze = make_shared<Analyze>(space, pair, fileName);
EOF

class=analyze
args="(space, pair, fileName)"
shrptr=yes
printFactory
appendPrint

cat << EOF >> $ccfile
#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_
EOF

