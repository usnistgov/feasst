#include <gtest/gtest.h>
#include <limits.h>
#include "functions.h"
#include <algorithm>
#include <fstream>
#include "mins.h"
#include <complex>
#include "base_math.h"

TEST(Functions, mySign) {
	EXPECT_EQ(-15, mySign(15,-1));
	EXPECT_EQ(-15.2, mySign(15.2,-10.0));
	EXPECT_EQ(25, mySign(25,23.1));
}

// Check myFill by filling with ones, and summing elements

TEST(Functions, myFill2d) {
  vector<vector<double> > x;
  const int dim1 = 111; const int dim2 = 222;
  x.resize(dim1, vector<double>(dim2));
  myFill(1., x);
  double sum = 0.;
  for (int i = 0; i < dim1; ++i) {
    for (int j = 0; j < dim2; ++j) {
      sum += x[i][j];
    }
  }
  EXPECT_EQ(dim1*dim2,sum);
}

TEST(Functions, myFill3d) {
  vector<vector<vector<double> > > x;
  const int dim1 = 111; const int dim2 = 22; const int dim3 = 33;
  x.resize(dim1, vector<vector<double> >(dim2, vector<double>(dim3)));
  myFill(1., x);
  double sum = 0.;
  for (int i = 0; i < dim1; ++i) {
    for (int j = 0; j < dim2; ++j) {
      for (int k = 0; k < dim3; ++k) {
        sum += x[i][j][k];
      }
    }
  }
  EXPECT_EQ(dim1*dim2*dim3,sum);
}

// Check myAdd by filling one vector with 1's, another with 2's, and checking for all 3's on the second vector, while the first vector remains unchanged

TEST(Functions, myAdd2d) {
  vector<vector<double> > x;
  vector<vector<double> > y;
  const int dim1 = 111; const int dim2 = 222;
  x.resize(dim1, vector<double>(dim2));
  y.resize(dim1, vector<double>(dim2));
  myFill(1., x);
  myFill(2., y);
  myAdd(x, y);
  for (int i = 0; i < dim1; ++i) {
    for (int j = 0; j < dim2; ++j) {
      EXPECT_EQ(3., y[i][j]);
      EXPECT_EQ(1., x[i][j]);
    }
  }
}

TEST(Functions, detTrSq) {
  vector<double> x(6);
	x[0] = 1;
	x[1] = 6;
	x[2] = 9;
	x[3] = 3;
	x[4] = 5;
	x[5] = 7;
	EXPECT_EQ(-16., det3DSym(x));
	EXPECT_EQ(16., tr3DSym(x));
  vector<double> y(6);
	sq3DSym(x, y);
	EXPECT_EQ(y[0], 35.);
	EXPECT_EQ(y[1], 94.);
	EXPECT_EQ(y[2], 155.);
	EXPECT_EQ(y[3], 56.);
	EXPECT_EQ(y[4], 71.);
	EXPECT_EQ(y[5], 120.);
}

// solve for eigen of symmetric matrix using cubicSolve
TEST(Functions, cubicSolve) {
  vector<double> x(6);
	x[0] = 1;
	x[1] = 6;
	x[2] = 9;
	x[3] = 3;
	x[4] = 5;
	x[5] = 7;
	const double tr = tr3DSym(x);
  vector<double> x2(6);
	sq3DSym(x, x2);
	const double tr2 = tr3DSym(x2);
	EXPECT_EQ(-16., det3DSym(x));
	EXPECT_EQ(16., tr);
	EXPECT_EQ(x2[0], 35.);
	EXPECT_EQ(x2[1], 94.);
	EXPECT_EQ(x2[2], 155.);
	EXPECT_EQ(x2[3], 56.);
	EXPECT_EQ(x2[4], 71.);
	EXPECT_EQ(x2[5], 120.);
	const double a = -tr;
	const double b = 0.5*(tr*tr - tr2);
	const double c = -det3DSym(x);
	//vector<double> y(3);
	double r1,r2,r3;
	cubicSolve(a, b, c, r1, r2, r3);
	EXPECT_NEAR(-1.4399, r2, 1e-4);
	EXPECT_NEAR(0.6623, r3, 1e-4);
	EXPECT_NEAR(16.7776, r1, 1e-4);
}

TEST(Functions, numLines) {
	EXPECT_EQ(266, numLines("../unittest/colMat1.txt"));
}

TEST(Functions, myMatMul3d) {
  vector<vector<double> > a;
  vector<double> x;
  vector<double> b;
  const int dim = 3;
  a.resize(dim, vector<double>(dim));
  x.resize(dim);
  b.resize(dim);
	a[0][0]=1;
	a[0][1]=2;
	a[0][2]=3;
	a[1][0]=4;
	a[1][1]=5;
	a[1][2]=6;
	a[2][0]=7;
	a[2][1]=8;
	a[2][2]=9;
	x[0]=1;
	x[1]=2;
	x[2]=3;
	myMatVecMul(a, x, b);
	EXPECT_NEAR(14, b[0], 1e-5);
	EXPECT_NEAR(32, b[1], 1e-5);
	EXPECT_NEAR(50, b[2], 1e-5);
	b.clear();
  b = myMatVecMul(a, x);
	EXPECT_NEAR(14, b[0], 1e-5);
	EXPECT_NEAR(32, b[1], 1e-5);
	EXPECT_NEAR(50, b[2], 1e-5);
}

TEST(Functions, myVecDotProd) {
  const int dimen = 3;
	vector<double> a(dimen);
  vector<double> b(dimen);
	fill(a.begin(), a.end(), 2.);
	fill(b.begin(), b.end(), 3.5);
	EXPECT_NEAR(21, myVecDotProd(a, b), 1e-10);
}

TEST(Functions, myVecCrosProd) {
  const int dimen = 3;
	vector<double> a(dimen);
  vector<double> b(dimen);
  vector<double> c(dimen);
	a[0] = 1;
	a[1] = 2;
	a[2] = 3;
	fill(b.begin(), b.end(), 3.5);
	myVecCrosProd(a, b, c);
	EXPECT_NEAR(-3.5, c[0], 1e-10);
	EXPECT_NEAR(7., c[1], 1e-10);
	EXPECT_NEAR(-3.5, c[2], 1e-10);
}

TEST(Functions, minimumImage) {

	const int dimen = 3;
	vector<double> xi(dimen);
	xi[0] = 31.4; xi[1] = 29.9; xi[2] = 11.97;
	vector<double> xj(dimen);
	xj[0] = -0.04; xj[1] = 13.72; xj[2] = 23.64;
	vector<double> xij(dimen);
	vector<vector<double> > boxMat;
  boxMat.resize(dimen, vector<double>(dimen, 0));
	boxMat[0][0] = 34.8579;
	boxMat[1][0] = 0;
	boxMat[2][0] = 0;
	boxMat[0][1] = 11.6193;
	boxMat[1][1] = 32.8644;
	boxMat[2][1] = 0;
	boxMat[0][2] = -11.6193;
	boxMat[1][2] = 16.4322;
	boxMat[2][2] = 28.4614;
	minimumImage(xi, xj, boxMat, xij);
	EXPECT_NEAR(8.20140, xij[0], 1e-4);
	EXPECT_NEAR(-0.2522, xij[1], 0.001);
	EXPECT_NEAR(16.7914, xij[2], 0.001);

// comments below as previously implemented in frame,
// and moved to functions
//  Frame TRJ(3);
//  TRJ.initFrameGRO("../test/1L2Y_prodp1_2.gro");
// 	TRJ.readGRO("../test/1L2Y_prodp1_2.gro");
//  const int dimen = 3;
//	vector<double> xij(dimen);
//	TRJ.minimumImage(TRJ.x[0], TRJ.x[1], xij);
//	EXPECT_NEAR(0.3, xij[0], 0.01);
//
//	TRJ.minimumImage(TRJ.x[0], TRJ.x[3199], xij);
//	EXPECT_NEAR(-3.4179, xij[0], 1e-4);
//	EXPECT_NEAR(16.18, xij[1], 0.001);
//	EXPECT_NEAR(-11.67, xij[2], 0.001);
//
//	TRJ.addGroup("atom", 0, 0, TRJ.natom() - 1);
//	TRJ.addGroupPDBComp("../test/1L2Y_prodp1_2.pdb", "test/1L2Y_prodp1_2ow.pdb");
//	TRJ.groupReduce();
//	double d;
//	for (int i = 0; i < TRJ.ngroupForType[1] - 1; ++i) {
//		for (int j = i + 1; j < TRJ.ngroupForType[1]; ++j) {
//			TRJ.minimumImage(TRJ.x[i], TRJ.x[j], xij);
//			d = pow(xij[0]*xij[0]+xij[1]*xij[1]+xij[2]*xij[2], 0.5);
//			
//		}
//	}
	
}

TEST(Functions, sphereGrid) {
	const int nTheta = 20;
	const int nPhi = 40;
	vector<vector<double> > x;
	x.resize(nTheta*nPhi, vector<double>(3, 0));
	sphereGrid(nTheta, nPhi, x);
	for (int i = 0; i < int(x.size()); ++i) {
		EXPECT_NEAR(1, myVecDotProd(x[i], x[i]), 1e-10);
	}
}

TEST(Functions, quat2rot) {
  vector<double> x;
  x.resize(4);
  x[0] = 0;
  x[1] = 0;
  x[2] = 0;
  x[3] = 1;
  vector<vector<double> > r;
  const int dim1 = 3; const int dim2 = 3;
  r.resize(dim1, vector<double>(dim2));
  quat2rot(x, r);
  EXPECT_EQ(r[0][0], 1);
  EXPECT_EQ(r[1][1], 1);
  EXPECT_EQ(r[2][2], 1);
  EXPECT_EQ(r[0][1], 0);
  EXPECT_EQ(r[0][2], 0);
  EXPECT_EQ(r[1][0], 0);
  EXPECT_EQ(r[1][2], 0);
  EXPECT_EQ(r[2][0], 0);
  EXPECT_EQ(r[2][1], 0);
  r.clear();
  r = quat2rot(x);
  EXPECT_EQ(r[0][0], 1);
  EXPECT_EQ(r[1][1], 1);
  EXPECT_EQ(r[2][2], 1);
  EXPECT_EQ(r[0][1], 0);
  EXPECT_EQ(r[0][2], 0);
  EXPECT_EQ(r[1][0], 0);
  EXPECT_EQ(r[1][2], 0);
  EXPECT_EQ(r[2][0], 0);
  EXPECT_EQ(r[2][1], 0);
}

TEST(Functions, myProd) {
  vector<double> x(4);
  x.at(0) = 1.5; x.at(1) = 2.3; x.at(2) = 3.5; x.at(3) = 50;
  EXPECT_EQ(1.5*2.3*3.5*50, myProd(x));
  vector<int> xi(4);
  x.at(0) = 1; x.at(1) = 2; x.at(2) = 3; x.at(3) = 50;
  EXPECT_EQ(300, myProd(x));
}

TEST(Functions, PI) {
  EXPECT_NEAR(3.14159265358979323846264338327950, PI, 1e-170);
}

TEST(Functions, myMatMul) {
  vector<vector<double> > a;
  vector<vector<double> > b;
  vector<vector<double> > c;
  const int dim = 3;
  a.resize(dim, vector<double>(dim));
  b.resize(dim, vector<double>(dim));
  c.resize(dim, vector<double>(dim));
	a[0][0]=1;
	a[0][1]=2;
	a[0][2]=3;
	a[1][0]=4;
	a[1][1]=5;
	a[1][2]=6;
	a[2][0]=7;
	a[2][1]=8;
	a[2][2]=9;
	b[0][0]=1;
	b[0][1]=0;
	b[0][2]=0;
	b[1][0]=0;
	b[1][1]=1;
	b[1][2]=0;
	b[2][0]=0;
	b[2][1]=0;
	b[2][2]=1;
	myMatMul(a, b, c);
	EXPECT_NEAR(c[0][0], 1, 1e-15);
	EXPECT_NEAR(c[0][1], 2, 1e-15);
	EXPECT_NEAR(c[0][2], 3, 1e-15);
	EXPECT_NEAR(c[1][0], 4, 1e-15);
	EXPECT_NEAR(c[1][1], 5, 1e-15);
	EXPECT_NEAR(c[1][2], 6, 1e-15);
	EXPECT_NEAR(c[2][0], 7, 1e-15);
	EXPECT_NEAR(c[2][1], 8, 1e-15);
	EXPECT_NEAR(c[2][2], 9, 1e-15);
  a.resize(2, vector<double>(3));
  b.resize(3, vector<double>(2));
  c.resize(2, vector<double>(2));
  a[0][0] = 1;
  a[0][1] = 0;
  a[0][2] = -2;
  a[1][0] = 0;
  a[1][1] = 3;
  a[1][2] = -1;
  b[0][0] = 0;
  b[0][1] = 3;
  b[1][0] = -2;
  b[1][1] = -1;
  b[2][0] = 0;
  b[2][1] = 4;
  myMatMul(a, b, c);
  EXPECT_NEAR(c[0][0], 0, 1e-18);
  EXPECT_NEAR(c[0][1], -5, 1e-18);
  EXPECT_NEAR(c[1][0], -6, 1e-18);
  EXPECT_NEAR(c[1][1], -7, 1e-18);
  c.clear();
  c = myMatMul(a, b);
  EXPECT_NEAR(c[0][0], 0, 1e-18);
  EXPECT_NEAR(c[0][1], -5, 1e-18);
  EXPECT_NEAR(c[1][0], -6, 1e-18);
  EXPECT_NEAR(c[1][1], -7, 1e-18);
}

TEST(Functions, vecSPCE) {
  vector<vector<double> > x = vecSPCE();
  EXPECT_EQ(0, x[0][0]);
  EXPECT_NEAR(-0.333313247568237, x[2][0], 1e-15);
  EXPECT_NEAR(0.942816142731718, x[2][1], 1e-15);
}

TEST(Functions, erfcc) {
  for (double x = -100; x <= 100; x += 0.05) {
    EXPECT_NEAR(erfc(x), erfcc(x), 1e-7);
  }
}

TEST(Functions, intersect) {
  vector<int> a, b;
  for (int i = 5; i < 30; i += 5) a.push_back(i);
  for (int i = 50; i > 9; i -= 10) b.push_back(i);
  vector<int> c(int(a.size()));
  //for (int i = 10; i < 51; i += 10) b.push_back(i);
  //for (int i = 0; i < int(a.size()); ++i) std::cout << a[i] << std::endl;
  //for (int i = 0; i < int(b.size()); ++i) std::cout << b[i] << std::endl;
  vector<int>::iterator it;
  std::sort (a.begin(), a.end());
  std::sort (b.begin(), b.end());
  //for (int i = 0; i < int(b.size()); ++i) std::cout << b[i] << std::endl;
  it=std::set_intersection (a.begin(), a.end(), b.begin(), b.end(), c.begin());
  //it=std::set_intersection (a.begin(), a.end(), b.begin(), b.end(), back_inserter(c));
  c.resize(it-c.begin());
  //std::cout << "The intersection has " << (c.size()) << " elements:\n";
  //for (it=c.begin(); it!=c.end(); ++it) {
  //  std::cout << ' ' << *it;
  //}
  //std::cout << '\n';
}

TEST(Functions, myVecAv) {
  vector<int> x;
  for (int i = 0; i < 6; ++i) {
    x.push_back(i);
  }
  EXPECT_NEAR(2.5, myVecAv(x), 1e-14);
}

TEST(Functions, fileOperations) {
  std::ofstream file("tmp/tmp.txt");
  file << "hi" << std::endl;
  EXPECT_TRUE(myFileExists("tmp/tmp.txt"));
  fileBackUp("tmp/tmp.txt");
  EXPECT_FALSE(myFileExists("tmp/tmp.txt"));
  EXPECT_TRUE(myFileExists("tmp/tmp.txt.bak"));
}

TEST(Functions, nWindow) {
  const int nMin = 0, nMax = 265, nThreads = 24, nOverlap = 2;
  const double nExp = 2;
  //const int nMin = 0, nMax = 265, nThreads = 24, nExp = 2, nOverlap = 2;
  //const int nMin = 0, nMax = 10000, nThreads = 8, nExp = 3, nOverlap = 2;
  vector<vector<int> > nRange = nWindow(nMin, nMax, nExp, nThreads, nOverlap);
  EXPECT_EQ(nThreads, int(nRange.size()));
  for (int i = 0; i < int(nRange.size()); ++i) EXPECT_EQ(2, int(nRange[i].size()));
  EXPECT_EQ(nMin, nRange[0][0]);
  EXPECT_EQ(nMax, nRange[nThreads-1][1]);
  EXPECT_EQ(95,  nRange[2][1]);
  EXPECT_EQ(143, nRange[7][0]);

  const double mMin = -0.005, mMax = 0.505;
  vector<vector<double> > win = nWindowGrowth(mMin, mMax, 0.0, 12, 0.01);
  EXPECT_EQ(12, int(win.size()));
  for (int i = 0; i < int(win.size()); ++i) EXPECT_EQ(2, int(win[i].size()));
//  for (int i = 0; i < int(win.size()); ++i) cout << win[i][0] << " " << win[i][1] << endl;
  EXPECT_EQ(mMin, win[0][0]);
  EXPECT_EQ(mMax, win[11][1]);
  EXPECT_NEAR(0.115, win[3][0], doubleTolerance);

  win = nWindowGrowth(mMin, mMax, 0.0, 12, 0.01, 3);
  EXPECT_EQ(12, int(win.size()));
  for (int i = 0; i < int(win.size()); ++i) EXPECT_EQ(2, int(win[i].size()));
//  for (int i = 0; i < int(win.size()); ++i) cout << win[i][0] << " " << win[i][1] << " " << win[i][1] - win[i][0] << endl;
  EXPECT_EQ(mMin, win[0][0]);
  EXPECT_EQ(mMax, win[11][1]);
  EXPECT_NEAR(0.095, win[3][0], doubleTolerance);
}

//TEST(Functions, pos2quat) {
//  myRanInitByDate();
//  for (int i = 0; i < 10; ++i) {
//    vector<vector<double> > x = myRanVecSPCE();
//    vector<double> q = pos2quat(x, vecSPCE());
//    vector<vector<double> > x2 = quat2pos(q, vecSPCE());
//    EXPECT_EQ(int(x.size()), int(x2.size()));
//    for (int i = 0; i < int(x.size()); ++i) {
//      EXPECT_EQ(int(x[i].size()), int(x2[i].size()));
//      for (int j = 0; j < int(x[i].size()); ++j) {
//        EXPECT_NEAR(x[i][j], x2[i][j], 1e-18);
//      }
//    }
//  }
//}

TEST(Functions, vecOnePatch) {
  vector<vector<double> > x = vecOnePatch();
  EXPECT_EQ(0, x[0][0]);
  EXPECT_EQ(0, x[0][1]);
  EXPECT_EQ(0, x[0][2]);
  EXPECT_EQ(0, x[1][1]);
  EXPECT_EQ(0, x[1][2]);
  EXPECT_NEAR(1, x[1][0], 1e-16);
}

double func(const double x) {
  return x*x;
}

TEST(Functions, mins) {
  Golden golden;
  golden.bracket(4, 5, func);
  double xmin = golden.minimize(func);
  EXPECT_NEAR(0, xmin, 1e-7);
}

TEST(Functions, myTrim) {
  const char* f = "/my/path/to/file";
  std::string fs = myTrim("/", f);
  EXPECT_EQ(0, fs.compare("file"));

  const char* f2 = "file";
  fs.assign(myTrim("/", f2));
  EXPECT_EQ(0, fs.compare("file"));

  const char* f3 = "/home/username/feasst/forcefield/cg7mabaniso.json";
  EXPECT_EQ("json", myTrim(".", f3));
}

TEST(Functions, findLocalMaxima) {
  vector<double> negxSqData;
  vector<double> poly3;
  for (int i = 0; i < 100; ++i) {
    double x = i - 50;
    negxSqData.push_back(-x*x);
//    cout << i << " " << negxSqData.back() << endl;
    poly3.push_back(x-10*x*x+x*x*x);
  }
  for (int tol = 1; tol < 6; ++tol) {
    vector<int> max = findLocalMaxima(negxSqData, tol);
    EXPECT_EQ(1, int(max.size()));
    EXPECT_EQ(50, max[0]);
    max = findLocalMaxima(poly3, 1);
    EXPECT_EQ(2, int(max.size()));
    EXPECT_EQ(50, max[0]);
    EXPECT_EQ(99, max[1]);
  }
}

TEST(Functions, mySq) {
  EXPECT_EQ(4, mySq(2));
  EXPECT_EQ(6.25, mySq(2.5));
}

TEST(Functions, myTrunc) {
  EXPECT_EQ(210, myTrunc(211, 10));
  EXPECT_EQ(210, myTrunc(219, 10));
  EXPECT_EQ(220, myTrunc(220, 10));
  EXPECT_EQ(220, myTrunc(229, 10));
}

TEST(Functions, fstos) {
  EXPECT_EQ(1, fstod("lnz", "../unittest/colMat6.txt"));
  EXPECT_EQ(9.53674e-07, fstod("lnf", "../unittest/colMat6.txt"));
  string strtmp = fstos("cellType", "../unittest/colMat6.txt");
  EXPECT_TRUE(strtmp.empty());
}

TEST(Functions, eigenANDrotateJacobi) {
  vector<vector<double> > matrix(3, vector<double>(3, 0.)), evectors = matrix;
  matrix[0][0] = 7;   matrix[1][0] = -2;  matrix[2][0] = 0;
  matrix[0][1] = -2;  matrix[1][1] = 6;   matrix[2][1] = -2;
  matrix[0][2] = 0;   matrix[1][2] = -2;  matrix[2][2] = 5;
  vector<double> evalues;
  jacobi(matrix, evalues, evectors);
  EXPECT_NEAR(evalues[0], 9, 1e-10);
  EXPECT_NEAR(evalues[1], 3, 1e-10);
  EXPECT_NEAR(evalues[2], 6, 1e-10);
  EXPECT_NEAR(evectors[0][0], 2/3., 1e-5);
  EXPECT_NEAR(evectors[1][0], -2/3., 1e-5);
  EXPECT_NEAR(evectors[2][0], 1/3., 1e-5);
  EXPECT_NEAR(evectors[0][1], 1/3., 1e-5);
  EXPECT_NEAR(evectors[1][1], 2/3., 1e-5);
  EXPECT_NEAR(evectors[2][1], 2/3., 1e-5);
  EXPECT_NEAR(evectors[0][2], -2/3., 1e-5);
  EXPECT_NEAR(evectors[1][2], -1/3., 1e-5);
  EXPECT_NEAR(evectors[2][2], 2/3., 1e-5);
}

TEST(Functions, sqDiff) {
  vector<double> x, y;
  for (int i = 0; i < 10; ++i) {
    x.push_back(i);
    y.push_back(i+1);
  }
  EXPECT_NEAR(10, sqDiff(x, y), 1e-19);
}

TEST(Functions, theta2rot) {
  const double theta = PI/8.;
  vector<vector<double> > r = theta2rot(theta);
  vector<double> unit(2);
  unit[0]=1;
  vector<double> ans = myMatVecMul(r, unit);
  EXPECT_NEAR(ans[0], cos(theta), 1e-19);
  EXPECT_NEAR(ans[1], sin(theta), 1e-19);

  double theta2 = theta+6*PI;
  theta2 += pbc2d(theta2);
  //theta2 -= int(theta2/2/PI)*2*PI;
  EXPECT_NEAR(theta, theta2, 1e-14);

}

TEST(Functions, orthogonalVec) {
  BaseMath bm;
  for (int i = 0; i < 100; ++i) {
    vector<double> xa = bm.ranUnitSphere(3);
    vector<double> xb = orthogonalVec(xa);
    EXPECT_NEAR(0, myVecDotProd(xa, xb), doubleTolerance);
  }
}

TEST(Functions, rotateVecByAxisAngle) {
  BaseMath bm;
//  myRanInitForRepro();
//  for (int i = 0; i < 1; ++i) {
//  //for (int i = 0; i < 100; ++i) {
//    // take a random unit vector and angle
//    vector<double> xa = bm.ranUnitSphere(3);
//    const double theta = PI*bm.uniformRanNum();
//    cout << "theta " << theta << " xa " << xa[0] << " " << xa[1] << " " << xa[2] << endl;
//
//    // rotate xa by theta about the z axis
//    vector<double> z(3); z[0]=0; z[1]=0; z[2]=1;
//    vector<double> xb = rotateVecByAxisAngle(xa, z, theta);
//    cout << "xb " << xb[0] << " " << xb[1] << " " << xb[2] << endl;
//
//    // use the law of cosines to compute the angle between A, origin, and B
//    const double xal = sqrt(myVecDotProd(xa, xa));
//    const double xbl = sqrt(myVecDotProd(xa, xa));
//    vector<double> xab(3); xab[0] = xa[0]-xb[0]; xab[1] = xa[1]-xb[1]; xab[2] = xa[2]-xb[2];
//    const double xabl = sqrt(myVecDotProd(xab, xab));
//    cout << "xal " << xal << " xbl " << xbl << " xab " << xabl << endl;
//
//    const double cosC = (xal*xal+xbl*xbl-xabl*xabl)/(2.*xal*xbl);
//
//    EXPECT_NEAR(cosC, cos(theta), doubleTolerance);
////    // project onto the xy-plane and compute the angle
////    const double projafac = myVecDotProd(xa, z);
////    const double projbfac = myVecDotProd(xb, z);
////    vector<double> xap = xa, xbp = xb;
////    for (int dim = 0; dim < 3; ++dim) {
////      xap[dim] = xa[dim] - projafac*z[dim];
////      xbp[dim] = xb[dim] - projbfac*z[dim];
////    }
////    EXPECT_NEAR(cos(180-theta), myVecDotProd(xap, xbp), 1e-19);
//  }

  vector<double> xa(3);
  xa[0] = 1; xa[1] = 0; xa[2] = 0;
  vector<double> z(3); z[0]=0; z[1]=0; z[2]=1;
  vector<double> xb = rotateVecByAxisAngle(xa, z, PI/2.);
  EXPECT_NEAR(0, xb[0], doubleTolerance);
  EXPECT_NEAR(1, xb[1], doubleTolerance);
  EXPECT_NEAR(0, xb[2], doubleTolerance);
}

TEST(Functions, quadraticEqReal) {
  double x1, x2;
  quadraticEqReal(1., -3., 2., x1, x2);
  EXPECT_NEAR(2, x1, doubleTolerance);
  EXPECT_NEAR(1, x2, doubleTolerance);
  quadraticEqReal(1., -11., 30., x1, x2);
  EXPECT_NEAR(6, x1, doubleTolerance);
  EXPECT_NEAR(5, x2, doubleTolerance);
  BaseMath bm;
  for (int i = 0; i < 100; ++i) {
    const double x1a = bm.uniformRanNum(), x2a = bm.uniformRanNum();
    quadraticEqReal(1., -(x1a+x2a), x1a*x2a, x1, x2);
    EXPECT_NEAR(x1a*x2a, x1*x2, doubleTolerance);
    EXPECT_NEAR(x1a+x2a, x1+x2, doubleTolerance);
  }
}

TEST(Functions, superQuadric) {
  const int n = 17;
  const double xdiff = 3;
  vector<double> axis(3); axis[2] = 1.;
  vector<vector<double> > sq = superQuadric3D(1., 3., n);
  vector<vector<double> > sq2 = superQuadric3D(1., 3., n);
  for (int i = int(sq.size()) - 1; i >= 0; --i) {
    if (sq[i][0] < 0) {
      sq.erase(sq.begin() + i);
    }
  }
  for (int i = int(sq2.size()) - 1; i >= 0; --i) {
    sq2[i] = rotateVecByAxisAngle(sq2[i], axis, PI/6);
    sq2[i][0] += xdiff;
    if (sq2[i][0] > xdiff) {
      sq2.erase(sq2.begin() + i);
    }
  }
//  for (int i = 0; i < int(sq.size()); ++i) {
//    for (int j = 0; j < int(sq[i].size()); ++j) {
//      cout << sq[i][j] << " ";
//    }
//    cout << endl;
//  }
//  for (int i = 0; i < int(sq2.size()); ++i) {
//    for (int j = 0; j < int(sq2[i].size()); ++j) {
//      cout << sq2[i][j] << " ";
//    }
//    cout << endl;
//  }
  double xmin;
  int x1min, x2min;
  minDistGrid(sq, sq2, xmin, x1min, x2min);
  EXPECT_NEAR(xdiff-2., xmin, 0.1);
  //EXPECT_NEAR(xdiff-2., xmin, sqrt(doubleTolerance));
  vector<double> xd(3);
  xd[0] = sq2[x2min][0] - sq[x1min][0];
  xd[1] = sq2[x2min][1] - sq[x1min][1];
  xd[2] = sq2[x2min][2] - sq[x1min][2];
  EXPECT_NEAR(xdiff-2, xd[0], 0.1);
  EXPECT_NEAR(0, xd[1], 0.02);
  EXPECT_NEAR(0, xd[2], doubleTolerance);
}

TEST(Functions, cartesian2spherical) {
  vector<double> x(3); x[0] = 1; x[1] = 0; x[2] = 0;
  vector<double> xS = cartesian2spherical(x);
  EXPECT_NEAR(xS[0], 1, doubleTolerance);
  EXPECT_NEAR(xS[1], PI/2, doubleTolerance);
  EXPECT_NEAR(xS[2], 0, doubleTolerance);

  vector<std::complex<double> > sphH = cart2sphereHarm6(x);
  const double q6 = sqrt(4*PI/13.* complexVec2norm(sphH));
  EXPECT_NEAR(1, q6, doubleTolerance);
}

TEST(Functions, Euler2RotMatANDRotMat2Euler) {
  BaseMath b;
  for (int trial = 0; trial < 100; ++trial) {
    // generate random euler angles
    vector<double> euler(3);
    euler[0] = (b.uniformRanNum()-0.5)*2*PI;
    euler[1] =  b.uniformRanNum()*PI;
    euler[2] = (b.uniformRanNum()-0.5)*2*PI;

    // test gimbal lock
    if (trial%2==0) euler[1] = 1e-15;

    vector<vector<double> > rotMat = Euler2RotMat(euler);
    vector<vector<double> > eulerConv = RotMat2Euler(rotMat);
    vector<vector<double> > rotMatConv = Euler2RotMat(eulerConv[0]);

    const double tol = 1e-10;
    //const double tol = sqrt(doubleTolerance);
    //EXPECT_NEAR(euler[0], eulerConv[0][0], tol);
    EXPECT_NEAR(euler[1], eulerConv[0][1], tol);
    //EXPECT_NEAR(euler[2], eulerConv[0][2], tol);
    for (int idim = 0; idim < 3; ++idim) {
      for (int jdim = 0; jdim < 3; ++jdim) {
        EXPECT_NEAR(rotMat[idim][jdim], rotMatConv[idim][jdim], tol);
      }
    }
  }
}

TEST(Functions, 3by3MatrixOperations) {
  vector<vector<double> > mat(3, vector<double>(3, 0));
  mat[0][0] = 1;
  mat[0][1] = 2;
  mat[0][2] = 0;
  mat[1][0] = -1;
  mat[1][1] = 1;
  mat[1][2] = 1;
  mat[2][0] = 1;
  mat[2][1] = 2;
  mat[2][2] = 3;
  EXPECT_NEAR(det3by3(mat), 9, doubleTolerance);
  vector<vector<double> > cof = cofactor3by3(mat);
  EXPECT_NEAR(1, cof[0][0], doubleTolerance);
  EXPECT_NEAR(4, cof[0][1], doubleTolerance);
  EXPECT_NEAR(-3, cof[0][2], doubleTolerance);
  EXPECT_NEAR(-6, cof[1][0], doubleTolerance);
  EXPECT_NEAR(3, cof[1][1], doubleTolerance);
  EXPECT_NEAR(0, cof[1][2], doubleTolerance);
  EXPECT_NEAR(2, cof[2][0], doubleTolerance);
  EXPECT_NEAR(-1, cof[2][1], doubleTolerance);
  EXPECT_NEAR(3, cof[2][2], doubleTolerance);
  const vector<vector<double> > inv = inv3by3(mat);
  EXPECT_NEAR(1./9, inv[0][0], doubleTolerance);
  EXPECT_NEAR(-6./9, inv[0][1], doubleTolerance);
  EXPECT_NEAR(2./9, inv[0][2], doubleTolerance);
  EXPECT_NEAR(4./9, inv[1][0], doubleTolerance);
  EXPECT_NEAR(3./9, inv[1][1], doubleTolerance);
  EXPECT_NEAR(-1./9, inv[1][2], doubleTolerance);
  EXPECT_NEAR(-3./9, inv[2][0], doubleTolerance);
  EXPECT_NEAR(0, inv[2][1], doubleTolerance);
  EXPECT_NEAR(3./9, inv[2][2], doubleTolerance);
}

TEST(Functions, myMinMaxElement) {
  vector<vector<vector<vector<vector<double> > > > > data(5, vector<vector<vector<vector<double> > > >(6, vector<vector<vector<double> > >(3, vector<vector<double> >(2, vector<double>(12, 1.2)))));
  data[1][2][1][0][4] = 5;
  data[4][0][2][1][9] = -1.;
//  EXPECT_NEAR(5, myMaxElement(data), doubleTolerance);
  EXPECT_NEAR(-1, myMinElement(data), doubleTolerance);
}

TEST(Functions, eulerPBC) {
  vector<double> euler(3);
  euler[0] = 60; euler[2] = -350; euler[1] = 0.2531235;
  vector<vector<double> > rotMat = Euler2RotMat(euler);
  eulerPBC(euler);
  EXPECT_NEAR(euler[0], 60-10*2*PI, doubleTolerance);
  EXPECT_NEAR(euler[1], 0.2531235, doubleTolerance);
  EXPECT_NEAR(euler[2], -350+112*PI, 100*doubleTolerance);

  vector<vector<double> > rotMat2 = Euler2RotMat(euler);
  for (int idim = 0; idim<3;++idim) {
    for (int jdim = 0; jdim<3;++jdim) {
      EXPECT_NEAR(rotMat[idim][jdim], rotMat2[idim][jdim], 100*doubleTolerance);
    }
  }
}

TEST(Functions, outerProd) {
  vector<double> u = {1, 2, 3};
  vector<double> v = {4, 5, 6};
  vector<vector<double> > w = outerProd(u, v);
  EXPECT_EQ(w[0][0], 4);
  EXPECT_EQ(w[0][1], 5);
  EXPECT_EQ(w[0][2], 6);
  EXPECT_EQ(w[1][0], 8);
  EXPECT_EQ(w[1][1], 10);
  EXPECT_EQ(w[1][2], 12);
  EXPECT_EQ(w[2][0], 12);
  EXPECT_EQ(w[2][1], 15);
  EXPECT_EQ(w[2][2], 18);
}

TEST(Functions, rotInversion) {
  vector<double> xref = {0, 0, 1.};
  vector<double> euler = {1.1, 2.2, 0.3};
  vector<vector<double> > rot = Euler2RotMat(euler);
  vector<double> x = myMatVecMul(rot, xref);
  vector<vector<double> > a1 = outerProd(x, xref);
  vector<vector<double> > a2 = outerProd(xref, xref);
  vector<vector<double> > a3 = transpose(a2);
  vector<vector<double> > rot2 = myMatMul(a1, a3);
  vector<double> x2 = myMatVecMul(rot2, xref);
  for (unsigned int dim = 0; dim < x.size(); ++dim) {
    EXPECT_NEAR(x[dim], x2[dim], doubleTolerance);
  }
}

