# This script makes the feasst.i, feasst.h and factories.cc files

from __future__ import print_function

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--source_dir", "-s", help="/path/to/feasst", default="err", type=str)
args = parser.parse_args()

assert(args.source_dir != "err")  # -s argument is required

factoryBaseClasses=["pair", "criteria", "analyze", "trial", "random"]
baseClasses=["accumulator", "analyze", "barrier", "base", "criteria", "custom_exception", "histogram", "mc", "pair", "random", "shape", "space", "table", "trial"]
nonClasses=["functions", "ui_abbreviated", "physical_constants"]

####################################################
print("Generating feasst.i and feasst.h")
####################################################

def strip_char(text, char="/", direction="back"):
    ''' use this function to remove path after glob'''
    assert(direction == "back")
    index = text.rfind(char)
    return text[(index+1):]
assert(strip_char("/home/user/hi") == "hi")

import glob
def globfiles(cls="", file_type=".h"):
    ''' use this function to return a list of files from src '''
    return sorted(glob.glob(args.source_dir+"/src/"+cls+"*"+file_type))

def printer(prepend, append, print_file, class_list=[""]):
    for cls in class_list:
        for full_header in globfiles(cls):
            header=strip_char(full_header)
            if header != "feasst.h":
                print(prepend+header+append, file=print_file)

import re
def my_grep(file_name, regex):
    lines = list()
    with open(file_name) as origin_file:
        for line in origin_file:
            line = re.findall(regex, line)
            #r'class (.*) : public (.*) {', line)
            if line:
                lines.append(line[0][0])
    return lines

def print_factory(print_file, base_class, append,
                  prepend="  } else if (typestr.compare(\"",
                  middle1="\") == 0) { ", middle2=""):
    for full_header in globfiles(base_class):
        cls_list = my_grep(file_name=full_header, regex=r'class (.*) : public (.*) {')
        if len(cls_list) != 0:
            for cls in cls_list:
                if middle1 == "":
                    print(prepend+cls+append, file=print_file)
                else:
                    print(prepend+cls+middle1+base_class+middle2+
                          cls+append, file=print_file)
    if middle1 != "":
        print("  } else { ASSERT(0, \"unrecognized "+base_class+"(\" << typestr << \") in factory\"); }\n\
  return "+base_class+";\n\
}\n", file=print_file)

with open(args.source_dir+"/src/feasst.h", "w") as headerfile:
    printer(prepend="#include \"", append="\"", print_file=headerfile)

with open("feasst.i", "w") as swigfile:
    print(
"/* NOTE: This file is created by tools/build/makeFactory.sh automatically\n\
   by cmake. Do not modify this file. */\n\
\n\
%module feasst\n\
\n\
%ignore CriteriaWLTMMC::lnPIrwsatwrap;\n\
%ignore CriteriaWLTMMC::lnPIrwnmxwrap;\n\
%ignore MC::boyleminwrap;\n\
\n\
%{", file=swigfile)

    printer(prepend="#include \"", append="\"", print_file=swigfile)

    print(
"#ifdef XDRFILE_H_\n\
  extern \"C\" {\n\
    #include \"xdrfile.h\"\n\
    #include \"xdrfile_xtc.h\"\n\
    #include \"xdrfile_trr.h\"\n\
  }\n\
#endif // XDRFILE_H_\n\
%}\n\
\n\
%include \"std_vector.i\"\n\
%include \"std_shared_ptr.i\"\n\
%include \"std_string.i\"\n\
%include \"std_iostream.i\"\n\
%template(IntVector) std::vector<int>;\n\
%template(DoubleVector) std::vector<double>;\n\
using namespace std;\n\
// using std::vector;\n\
\n\
%pythonnondynamic;\n\
\n\
%shared_ptr(std::ifstream);\n\
namespace std{\n\
  class ifstream {\n\
  };\n\
}\n\
\n\
%include \"std_map.i\"\n\
\n\
%template(args) std::map<std::string, std::string>;\n\
\n\
%shared_ptr(Base);", file=swigfile)

    print_factory(print_file=swigfile, base_class="", prepend="%shared_ptr(",
                  append=");", middle1="")

    printer(prepend="%include ", append="", print_file=swigfile)

####################################################
print("Generating factories.cc")
####################################################

factories_file_name = args.source_dir+"/src/factories.cc"
with open(factories_file_name, "w") as factory_file:
    print(
"/* Factory method for Pair, Criteria, Trial and Analyze.\n\
   This file is created by tools/build/makeFactory.sh automatically\n\
   by cmake. Do not modify this file. */", file=factory_file)

    printer(prepend="#include \"", append="\"", print_file=factory_file,
            class_list=factoryBaseClasses)

    print(
"#ifdef FEASST_NAMESPACE_\n\
namespace feasst {\n\
#endif  // FEASST_NAMESPACE_\n\
\n\
Pair* makePair(Space* space, const char* fileName) {\n\
  ASSERT(fileExists(fileName),\n\
    \"restart file(\" << fileName << \") doesn't exist\");\n\
  string typestr = fstos(\"className\", fileName);\n\
  Pair* pair_ = NULL;\n\
  if (1 == 0) { // filler for a loop of \"} else if {\"", file=factory_file)

    print_factory(print_file=factory_file,
                  base_class="pair_",
                  middle2=" = new ",
                  append="(space, fileName);")

    print(
"Criteria* makeCriteria(const char* fileName) {\n\
  ASSERT(fileExists(fileName),\n\
    \"restart file(\" << fileName << \") doesn't exist\");\n\
  string typestr = fstos(\"className\", fileName);\n\
  Criteria* criteria_ = NULL;\n\
  if (typestr.compare(\"Criteria\") == 0) {\n\
    criteria_ = new CriteriaMetropolis(fileName);", file=factory_file)

    print_factory(print_file=factory_file,
                  base_class="criteria_",
                  middle2=" = new ",
                  append="(fileName);")

    print(
"shared_ptr<Trial> makeTrial(Pair* pair, Criteria* criteria,\n\
  const char* fileName) {\n\
  ASSERT(fileExists(fileName),\n\
         \"restart file(\" << fileName << \") doesn't exist\");\n\
  string typestr = fstos(\"className\", fileName);\n\
  shared_ptr<Trial> trial_;\n\
  if (1 == 0) { // filler for a loop of \"} else if {\"", file=factory_file)

    print_factory(print_file=factory_file,
                  base_class="trial_",
                  middle2=" = make_shared<",
                  append=">(fileName, pair, criteria);")

    print(
"shared_ptr<Analyze> makeAnalyze(Pair* pair,\n\
  const char* fileName) {\n\
  ASSERT(fileExists(fileName),\n\
         \"restart file(\" << fileName << \") doesn't exist\");\n\
  string typestr = fstos(\"className\", fileName);\n\
  shared_ptr<Analyze> analyze_;\n\
  if (typestr.compare(\"Analyze\") == 0) {\n\
    analyze_ = make_shared<Analyze>(pair, fileName);", file=factory_file)

    print_factory(print_file=factory_file,
                  base_class="analyze_",
                  middle2=" = make_shared<",
                  append=">(pair, fileName);")

    print(
"#ifdef FEASST_NAMESPACE_\n\
}  // namespace feasst\n\
#endif  // FEASST_NAMESPACE_", file=factory_file)



