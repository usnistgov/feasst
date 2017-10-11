/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include <gtest/gtest.h>
#include "table.h"
#include "functions.h"

#ifdef FEASST_NAMESPACE_
using namespace feasst;
#endif  // FEASST_NAMESPACE_

TEST(Table, roundSquare) {
  Table vt("../unittest/table/roundSquare/vtrg0.04nt11nz11nd6");
  Table rm("../unittest/table/roundSquare/rmrg0.04nt11nz11nd6");
  Table rc("../unittest/table/roundSquare/rcrg0.04nt11nz11nd6");

  EXPECT_NEAR(1, rm.interpolate(0.5, PI/4, PI/4), 3e-8);
  EXPECT_NEAR(1.08, rc.interpolate(0.5, 0,0), 3e-8);
  EXPECT_NEAR(-4.28831750, vt.interpolate(0.5, 0,0,0), 1e-10);
  EXPECT_NEAR(-4.28831750, vt.interpolate(0.5, 0,PI/4, PI/4), 1e-10);
  EXPECT_NEAR((-3.65815655+-4.28831750)/2., vt.interpolate(0.5, 0.05,PI/4, PI/4), 3e-8);

  EXPECT_NEAR(sqrt(PI/4), rm.interpolate(0., PI/4, PI/4), 1.2e-8);
  EXPECT_NEAR(sqrt(PI/2)+0.08, rc.interpolate(0., 0,0), 4e-8);
  EXPECT_NEAR(-1.13581623, vt.interpolate(0., 0,0,0), 1e-10);
  EXPECT_NEAR(-15.10313215, vt.interpolate(0., 0,PI/4, PI/4), 1e-10);
  EXPECT_NEAR((-15.10313215+-13.56564743)/2, vt.interpolate(0., 0.05,PI/4, PI/4), 3e-2);

  EXPECT_NEAR((0.88622693+0.89105682)/2, rm.interpolate(0.05, PI/4, PI/4), 1.2e-8);

  vt.compute_min_compress1d(0);
  EXPECT_NEAR(-15.10313215, vt.min(0), DTOL);
  EXPECT_NEAR(-4.28831750, vt.min(0.5), DTOL);
  EXPECT_NEAR((-13.13564708+-15.10313215)/2, vt.min(0.05), DTOL);
}

TEST(Table, superball) {
  Table rm("../unittest/table/superball/abcde0.5/k4nz3/tHard");

  EXPECT_NEAR(1., rm.interpolate(0,0,0,0,0), 2e-5);
  //EXPECT_NEAR(1.73204, rm.interpolate(0,0,4.18879,4.18879,1.5708), 2e-5);
  //EXPECT_NEAR(1.46411, rm.interpolate(0,3.14159,6.28319,4.18879,1.5708), 2e-5);
  //EXPECT_NEAR(0.5*(2+1.46411), rm.interpolate(0,3.14159,6.28319,4.18879,(3.14159+1.5708)/2), 2e-5);
}

TEST(Table, solidRevDisck6nz1) {
  Table rm("../unittest/table/solidRev/disc_k6nz1/tHard");
#ifdef HDF5_
  rm.printHDF5("tmp/table.hdf5");
#endif  //HDF5_
//  Table rm2("../unittest/table/solidRev/cylinder_ar1k150nz4/tHard");
//  rm2.printHDF5("tmp/tablecyl1.hdf5");
//  Table dep("../unittest/table/solidRev/cylinder_ar1k150nz4/tDep");
//  dep.printHDF5("tmp/tablecyl1dep.hdf5");
//  exit(0);
}

TEST(Table, lj) {
  Table lj("../unittest/lj/tabi0j0n100");
  //Table lj("../test/lj/table/tabi0j0n10");
  const double dr = 5e-2;
  //for (double r2 = 0.8*0.8; r2 <= 3*3.001; r2 += 0.836) {
    //cout << "r " << sqrt(r2) << " r2 " << r2 << " " << lj.interpolate(r2) << endl;
  //for (double r = 0.8; r <= 3.001; r += 0.22) {
  for (double r = 0.8; r <= 3+dr/2.; r += dr) {
    lj.interpolate(r*r);
    //const double pe = lj.interpolate(r*r);
    //cout << "r " << r << " " << lj.interpolate(r*r) << endl;
  }
}
