#include <gtest/gtest.h>
#include "pair_tabular.h"
#include "analyze_orient.h"

using namespace feasst;

TEST(AnalyzeOrient, orient) {
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(20,dim);
  s.addMolInit("../forcefield/data.onePatch1000z");
  s.initEuler(1);
  stringstream tableFile;
  tableFile << s.install_dir() << "/unittest/table/solidRev/spherek2nz4/rst";
  //tableFile << s.install_dir() << "/unittest/table/solidRev/cylinder_ar6k150nz4/rst";
  PairTabular p(&s, tableFile.str().c_str());
  p.initLMPData("../forcefield/data.onePatch1000z");
  p.initCuts();
  std::ifstream file("../unittest/solidRev/cylinder_ar6k150nz4/moviep0pr.xyz");
  p.readxyzeuler(file);

  AnalyzeOrient an(&s, &p);
  an.zbin=0.5;
  an.update();
  an.initFileName("tmp/odf.txt");
  an.print();
}


