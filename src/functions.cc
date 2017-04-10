/**
 * \file
 *
 * \brief external function library
 *
 * Implementation of utility functions for multidimensional vectors
 */

#include "functions.h"
#include <fstream>
#include <algorithm>
#include <math.h>
#include <complex>
#include <sys/stat.h>
#include <limits>
#include <signal.h>

/**
 * sum and difference operator
 */
double op_sum(double i, double j) { return i + j; }

/**
 * function to return the magnitude of the first argument, with the sign of the second argument
 */
double mySign(const double a,				//!< magnitude of return
	   					const double b				//!< sign of return
	) {

  return (b >= 0.0 ? fabs(a) : -fabs(a));
}

/**
 * Add two vectors of double vectors
 * The first argument isn't modified, second argument has the first added to it
 */
void myAdd(vector<vector<double> > &x,	//!< vector to add, unchanged
	   vector<vector<double> > &y		//!< vector added to, returns sum
	) {

  // check if lengths are equal
  if (x.size() != y.size()) myOutF("error", std::ostringstream().flush() << "x size(" << x.size() << ") != y size (" << y.size() << ").!");

  for (vector<vector<double> >::iterator xiter = x.begin(), yiter = y.begin();
       xiter != x.end() || yiter != y.end(); ++xiter,++yiter) {
    std::transform((*yiter).begin(), (*yiter).end(), (*xiter).begin(), (*yiter).begin(), op_sum);
  }

  return;
}

/**
 * Fill a vector of double vectors with a double
 */
void myFill(const double y, 	//!< value to fill vector
	    vector<vector<double> > &x	//!< vector to be filled
	) {

  for (vector<vector<double> >::iterator iter = x.begin();
       iter != x.end(); ++iter) {
    std::fill((*iter).begin(), (*iter).end(), y);
  }
}

/**
 * Fill a vector of int vectors with an int
 */
void myFill(const int y, 	//!< value to fill vector
	    vector<vector<int> > &x	//!< vector to be filled
	) {

  for (vector<vector<int> >::iterator iter = x.begin();
       iter != x.end(); ++iter) {
    std::fill((*iter).begin(), (*iter).end(), y);
  }
}

/**
 * Return determinant of symmetric 3d matrix given by 6d vector
 */
double det3DSym(const vector<double> &x	//!< matrix as 6d vec
	) {
	return x[0] * (x[1] * x[2] - x[5] * x[5]) -
				 x[3] * (x[3] * x[2] - x[5] * x[4]) +
				 x[4] * (x[3] * x[5] - x[1] * x[4]);
}

/**
 * Return trace of symmetric 3d matrix given by 6d vector
 */
double tr3DSym(const vector<double> &x	//!< matrix as 6d vec
	) {
	return x[0] + x[1] + x[2];
}

/**
 * Return square of 3d symmetric matrix given by 6d vector
 */
void sq3DSym(vector<double> x,	//!< matrix as 6d vec
               vector<double> &y	//!< matrix as 6d vec
	) {
	
	y[0] = x[0]*x[0] + x[3]*x[3] + x[4]*x[4];
	y[3] = x[0]*x[3] + x[3]*x[1] + x[5]*x[4];
	y[4] = x[0]*x[4] + x[5]*x[3] + x[2]*x[4];
	y[1] = x[3]*x[3] + x[1]*x[1] + x[5]*x[5];
	y[5] = x[3]*x[4] + x[1]*x[5] + x[2]*x[5];
	y[2] = x[4]*x[4] + x[5]*x[5] + x[2]*x[2];

}

/**
 * Return real roots of cubic equation
 */
void cubicSolve(const double b,		//!< cubic equation x^3 + bx^2 + cx + d
                const double c,		//!< cubic equation x^3 + bx^2 + cx + d
                const double d,		//!< cubic equation x^3 + bx^2 + cx + d
                double &r1,	//!< roots to cubic equation
                double &r2,	//!< roots to cubic equation
                double &r3 	//!< roots to cubic equation
	) {
	const double q = (3. * c - b * b) / 9.;
	const double r = (9. * b * c - 27. * d -
	                  2. * b * b * b) / 54.;
	const double v = q * q * q + r * r;
	std::complex<double> s, t;
	if ( v < 0.) {
	  s = pow(std::complex<double>(r, sqrt(-v)), 1./3.);
	  t = pow(std::complex<double>(r,-sqrt(-v)), 1./3.);
	} else {
	  s = pow(std::complex<double>(r + sqrt(v), 0.), 1./3.);
	  t = pow(std::complex<double>(r - sqrt(v), 0.), 1./3.);
	}
	r1 = -b/3. + s.real() + t.real();
	r2 = -b/3. - 0.5*(s.real() + t.real()) - 0.5*(s.imag() - t.imag())*sqrt(3);
	r3 = -b/3. - 0.5*(s.real() + t.real()) + 0.5*(s.imag() - t.imag())*sqrt(3);
}

/**
 * Return number of lines in file
 */
int numLines(const std::string fileName 			//!< file name
	) {

	// open file
  std::ifstream file(fileName.c_str());
  if (!file.good()) myOutF("error", std::ostringstream().flush() << "cannot open file " << fileName << " in functions.cc");

	int n = 0;
	std::string line;
	while(std::getline(file, line)) ++n;
	return n;
}

/**
 * Return product of matrix and vector
 *  A[][] x[] = b[]
 */
void myMatVecMul(vector<vector<double> > &a,		//!< coefficient matrix
								 vector<double> &x, 									//!< vector
								 vector<double> &b										//!< resulting vector
	) {

	// check dimensions
  if (int(a[0].size()) != int(x.size())) myOutF("error", std::ostringstream().flush() << "myMatMul requires columns of matrix a(" << int(a[0].size()) << ") to equal rows of vector x(" << int(x.size()));
  if (int(a.size()) != int(b.size())) myOutF("error", std::ostringstream().flush() << "myMatMul requires rows of matrix a(" << int(a.size()) << ") to equal rows of vector b(" << int(b.size()));

	//zero b vector
	fill(b.begin(), b.end(), 0.);

	// multiply
	for (int j = 0; j < int(a.size()); ++j) {
		for (int i = 0; i < int(x.size()); ++i) {
			b[j] += a[j][i] * x[i];
		}
	}
}

/**
 * Return product of matrix and vector
 *  A[][] x[] = b[]
 */
vector<double> myMatVecMul(vector<vector<double> > &a,		//!< coefficient matrix
								 vector<double> &x 									//!< vector
	) {

  vector<double> b(int(a.size()));

	// check dimensions
  if (int(a[0].size()) != int(x.size())) myOutF("error", std::ostringstream().flush() << "myMatMul requires columns of matrix a(" << int(a[0].size()) << ") to equal rows of vector x(" << int(x.size()));

	// multiply
	for (int j = 0; j < int(a.size()); ++j) {
		for (int i = 0; i < int(x.size()); ++i) {
			b[j] += a[j][i] * x[i];
		}
	}

  return b;
}

/**
 * Return vector dot product
 *  a.b = c
 */
double myVecDotProd(vector<double> &a,
								  vector<double> &b
	) {

	// check dimensions
  if (int(a.size()) != int(b.size())) myOutF("error", std::ostringstream().flush() << "myVecDotProd requires vectors of equal size, however, a(" << int(a.size()) << ") and b(" << int(b.size()));
	
	double c = 0;

	for (int i = 0; i < int(a.size()); ++i) {
		c += a[i]*b[i];
	}

	return c;
}

/*
 * function to normalize vector to unit size
 */
void myVecNormalize(vector<double> &x) {
  const double len = sqrt(myVecDotProd(x, x));
  for (int dim = 0; dim < int(x.size()); ++dim) {
    x[dim] /= len;
  }
}

/**
 * Return vector cross product
 *  a*b = c
 */
void myVecCrosProd(vector<double> &a,
							     vector<double> &b,
								   vector<double> &c
	) {

	// check dimensions
  if ( (int(a.size()) != int(b.size())) || (int(a.size()) != 3) ) myOutF("error", std::ostringstream().flush() << "myVecDotProd requires vectors of 3 dimensions, however, a(" << int(a.size()) << ") and b(" << int(b.size()));

	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

/**
 * return minimum image separation vector, xij
 *  xij = xi - xj
 */
void minimumImage(const vector<double> &xi,	//!< position vector of particle i
							    			 const vector<double> &xj, //!< position vector of particle j
							    			 const vector<vector<double> > &boxMat, //!< matrix defining boundary conditions
								  			 vector<double> &xij	//!< separation vector between i,j
	) {

	const int dimen = int(xi.size());
	double d2min = 0;
	double dx;
	vector<double> xijmin(dimen);
	vector<vector<double> > bm = boxMat;
	for (int dim = 0; dim < dimen; ++dim) {
		xijmin[dim] = dx = xi[dim] - xj[dim];
		d2min += dx*dx;
	}

	// brute force (slow) implementation
	// loop through all possible images and select the minimum
	double d2;
	vector<double> b(dimen);
	vector<double> n(dimen);
	for (int ix = -1; ix <= 1; ++ix) {
  for (int iy = -1; iy <= 1; ++iy) {
  for (int iz = -1; iz <= 1; ++iz) {
	  if ( (ix * ix + iy * iy + iz * iz) != 0) {
			
			// define shift vector n for displacements to periodic image box
			b[0]=ix; b[1]=iy; b[2]=iz;
			myMatVecMul(bm, b, n);

			// calculate separation vector to periodic image
			d2 = 0;
			for (int dim = 0; dim < dimen; ++dim) {
				dx = xi[dim] - xj[dim] + n[dim];
				d2 += dx*dx;
			}

			if (d2 < d2min) {
//				cout << "d2 " << d2 << " d2min " << d2min << endl;
//				cout << n[0] << " " << n[1] << " " << n[2] << endl;
//				cout <<  " xi " << xi[0] << " " << xi[1]  << " " << xi[2] << " xj " << xj[0] << " " << xj[1] << " " << xj[2] << endl;
//				cout <<  " xij " << xij[0] << " " << xij[1]  << " " << xij[2] << endl;
//				cout <<  " xijmin " << xijmin[0] << " " << xijmin[1]  << " " << xijmin[2] << endl;
				for (int dim = 0; dim < dimen; ++dim) xijmin[dim] = xi[dim] - xj[dim] + n[dim];
				d2min = d2;
			}
		}
	}}}
	xij = xijmin;

// below are two broken implementations as taken frome lammps source code

//	double xijdim;
//	double dx = xi[0]-xj[0];
//	double dy = xi[1]-xj[1];
//	double dz = xi[2]-xj[2];
//	if (fabs(dz) > 0.5 * boxMat[2][2]) {
//		dz -= mySign(boxMat[2][2], dz);
//		dy -= mySign(boxMat[2][1], dz);
//		dx -= mySign(boxMat[2][0], dz);
//	}
//	if (fabs(dy) > 0.5 * boxMat[1][1]) {
//		dy -= mySign(boxMat[1][1], dy);
//		dx -= mySign(boxMat[1][0], dy);
//	}
//	if (fabs(dx) > 0.5 * boxMat[0][0]) {
//		dx -= mySign(boxMat[0][0], dx);
//	}
//	xij[0]=dx;	
//	xij[1]=dy;	
//	xij[2]=dz;	

//	double xijdim;
//	for (int dim = 0; dim < int(xi.size()); ++dim) {
//		xijdim = xi[dim] - xj[dim];
//		if (fabs(xijdim) > 0.5 * boxMat[dim][dim]) {
//			for (int idim = 0; idim < int(xi.size()); ++idim) {
//				xijdim -= mySign(boxMat[dim][idim], xijdim);
//			}
//		}
//		xij[dim]=xijdim;
//	}
}

/**
 *  Generate a grid of points on a sphere, stored as cartesian coordinates in the vector x
 *   the grid is defined as the 2-d surface of nTheta equaly spaced values of theta
 *   and nPhi equally spaced values of phi
 */
void sphereGrid(const int nTheta, 	//!< number of points in theta
							  const int nPhi, 		//!< number of points in phi
							  vector<vector<double> > &x //!< output grid positions
	) {

	const double dtheta = PI/(nTheta - 1);
	const double dphi = 2.*PI/(nPhi - 1);
	int i = 0;
	for (int itheta = 0; itheta < nTheta; ++itheta) {
	  const double theta = dtheta*itheta;
		for (int iphi = 0; iphi < nPhi; ++iphi) {
			const double phi = dphi*iphi;
			
			// convert to cartesian coordinates
			x[i][2] = cos(theta);
			const double sintheta = sin(theta);
			x[i][0] = sintheta*cos(phi);
			x[i][1] = sintheta*sin(phi);

			++i;
		}	
	}
}

/**
 *  initialize random number generator based on date
 */
void myRanInitByDate() {
	const int t = time(NULL);
  srand ( t );
  myOutF("note", std::ostringstream().flush() << "time(seed): " << t);
}

/**
 *  initialize random number generator to seed value for reproducibility
 */
void myRanInitForRepro(const int seed) {
	srand ( seed );
  myOutF("warning", std::ostringstream().flush() << "Initializing random number generator for reproductionseed(" << seed << ")");
}

/**
 *  initialize random number generator to same value for reproducibility
 */
void myRanInitForRepro() {
	srand ( 1346867550 );
  myOutF("warning", std::ostringstream().flush() << "Initializing random number generator for reproduction seed(" << 1346867550 << ")");
}

/**
 * generate a rotation matrix from quaternions
 *  FRANZ J. VESELY, JOURNAL OF COMPUTATIONAL PHYSICS 47, 291-296 (1982)
 */
void quat2rot(vector<double> &q,	//!< quaterion vector of dimension D+1
	   vector<vector<double> > &r		//!< rotation matrix of size DxD
	) {

  // check dimensions
  if (q.size() != 4) myOutF("error", std::ostringstream().flush() << "3must be three dimensional, such that quaterion is a 4d vector\n");
  if (r.size() != 3) myOutF("error", std::ostringstream().flush() << "r must be three dimensional\n");
  if (r[0].size() != 3) myOutF("error", std::ostringstream().flush() << "r2 must be three dimensional\n");

  r[0][0] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
  r[1][0] = 2*(q[0]*q[1] + q[2]*q[3]);
  r[2][0] = 2*(q[2]*q[0] - q[1]*q[3]);
  r[0][1] = 2*(q[0]*q[1] - q[2]*q[3]);
  r[1][1] = q[1]*q[1] - q[2]*q[2] - q[0]*q[0] + q[3]*q[3];
  r[2][1] = 2*(q[1]*q[2] + q[0]*q[3]);
  r[0][2] = 2*(q[2]*q[0] + q[1]*q[3]);
  r[1][2] = 2*(q[1]*q[2] - q[0]*q[3]);
  r[2][2] = q[2]*q[2] - q[0]*q[0] - q[1]*q[1] + q[3]*q[3];
}

/**
 * generate a rotation matrix from quaternions
 *  FRANZ J. VESELY, JOURNAL OF COMPUTATIONAL PHYSICS 47, 291-296 (1982)
 */
vector<vector<double> > quat2rot(vector<double> q	//!< quaterion vector of dimension D+1
	) {

  vector<vector<double> > r(3, vector<double>(3));

  // check dimensions
//  if (q.size() != 4) {
//    raise(SIGSEGV);
//    myOutF("error", std::ostringstream().flush() << q[50]/0.<<"1must be three dimensional, such that quaterion is a 4d vector\n");
//  }
  ASSERT(int(q.size()) == 4, "1quaterion is a 4d vector in 3d");
  if (r.size() != 3) myOutF("error", std::ostringstream().flush() << "2must be three dimensional, such that quaterion is a 4d vector\n");
  if (r[0].size() != 3) myOutF("error", std::ostringstream().flush() << "3must be three dimensional, such that quaterion is a 4d vector\n");

  r[0][0] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
  r[1][0] = 2*(q[0]*q[1] + q[2]*q[3]);
  r[2][0] = 2*(q[2]*q[0] - q[1]*q[3]);
  r[0][1] = 2*(q[0]*q[1] - q[2]*q[3]);
  r[1][1] = q[1]*q[1] - q[2]*q[2] - q[0]*q[0] + q[3]*q[3];
  r[2][1] = 2*(q[1]*q[2] + q[0]*q[3]);
  r[0][2] = 2*(q[2]*q[0] + q[1]*q[3]);
  r[1][2] = 2*(q[1]*q[2] - q[0]*q[3]);
  r[2][2] = q[2]*q[2] - q[0]*q[0] - q[1]*q[1] + q[3]*q[3];

  return r;
}

/**
 * generate a 2d rotation matrix from angle
 */
vector<vector<double> > theta2rot(double theta) {
  vector<vector<double> > r(2, vector<double>(2));
  r[0][0] = cos(theta);
  r[0][1] = -sin(theta);
  r[1][0] = -r[0][1];
  r[1][1] = r[0][0];
  return r;
}

/**
 *  function to return vector x with atomic positions of a reference SPC/E water molecule
 *   vector of atoms, oxygen first
 */
vector<vector<double> > vecSPCE() {

  vector<vector<double> > x(3, vector<double>(3));
  const double doh = 1.;
  const double theta = 109.47;
  x[0][0] = 0.; x[0][1] = 0.; x[0][2] = 0;
  x[1][0] = doh; x[1][1] = 0.; x[1][2] = 0;
  x[2][0] = -1*doh*cos((180-theta)/180*PI);
  x[2][1] = doh*sin((180-theta)/180*PI);
  x[2][2] = 0;
  return x;
}

/**
  * function to return vector x with atomic positions of a reference one-patch molecule
 * */
vector<vector<double> > vecOnePatch() {
  vector<vector<double> > x(2, vector<double>(3, 0.));
  x[1][0] = 1.;
  return x;
}

/**
 * Return product of two matrices
 *  A[][] b[][] = c[][]
 */
void myMatMul(vector<vector<double> > &a,		//!< coefficient matrix
						  vector<vector<double> > &b,
							vector<vector<double> > &c
	) {

	// check dimensions
  if (int(a[0].size()) != int(b.size())) myOutF("error", std::ostringstream().flush() << "myMatMul requires columns of matrix a(" << int(a[0].size()) << ") to equal rows of vector x(" << int(b.size()));
  if (int(a.size()) != int(c.size())) myOutF("error", std::ostringstream().flush() << "myMatMul requires rows of matrix a(" << int(a.size()) << ") to equal rows of vector b(" << int(c.size()));

	//zero c
	myFill(0, c);

	// multiply
	for (int i = 0; i < int(a.size()); ++i) {
		for (int j = 0; j < int(b[0].size()); ++j) {
		  for (int k = 0; k < int(b.size()); ++k) {
		    c[i][j] += a[i][k]*b[k][j];	
		  }
    }
	}
}

/**
 * Return product of two matrices
 *  A[][] b[][] = c[][]
 */
vector<vector<double> > myMatMul(vector<vector<double> > &a,		//!< coefficient matrix
						  vector<vector<double> > &b
	) {

  vector<vector<double> > c(int(a.size()), vector<double>(int(b[0].size())));

	// check dimensions
  if (int(a[0].size()) != int(b.size())) myOutF("error", std::ostringstream().flush() << "myMatMul requires columns of matrix a(" << int(a[0].size()) << ") to equal rows of vector x(" << int(b.size()));
  if (int(a.size()) != int(c.size())) myOutF("error", std::ostringstream().flush() << "myMatMul requires rows of matrix a(" << int(a.size()) << ") to equal rows of vector b(" << int(c.size()));

	// multiply
	for (int i = 0; i < int(a.size()); ++i) {
		for (int j = 0; j < int(b[0].size()); ++j) {
		  for (int k = 0; k < int(b.size()); ++k) {
		    c[i][j] += a[i][k]*b[k][j];	
		  }
    }
	}

  return c;
}

/**
 * function to return complimentary error function, from Numerical Recipes 6.2 pg 221
 *  should be correct to 1.2e-7
 *  turns out it is slower than erfc in the math library, and less accurate
 */
double erfcc(const double x) {

  double t,z,ans;
  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
  t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
  t*(-0.82215223+t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
}

/**
 * function to return volume of spherical shell defined by rabove and rebelow in dim dimensions
 */
double volShell(const double rabove,     //!< upper limit of shell
                const double rbelow,     //!< lower limit of shell
                const int dim            //!< spatial dimensions
  ) {

  if (dim != 3) myOutF("error", std::ostringstream().flush() << "only implemented for 3 dimensions (dim=" << dim << ") using quaterions");
  return 4./3.*PI*(rabove*rabove*rabove - rbelow*rbelow*rbelow);
}

/**
 * function to return average of vector
 */
double myVecAv(const vector<int> &x
  ) {
  return double(std::accumulate(x.begin(), x.end(), 0)) / double(x.size());
  //return std::accumulate(x.begin(), x.end(), 0.) / int(x.size());
}

/**
 * function to create back-up
 */
bool myFileExists(const char* fileName    //!< file name in c_str
  ) {

  struct stat buf;
  if (stat(fileName, &buf) != -1) return true;
  return false;
}

/**
 * function to create back-up file
 */
void fileBackUp(const char* fileName    //!< file name in c_str
  ) {

  if (myFileExists(fileName)) {
    std::ostringstream f;
    f << fileName << ".bak";
    rename(fileName, f.str().c_str());
  }
}

/**
 * test the speed of arrays
 */
void arraySpeedTest() {
  vector<double> x(30000,1.);
  vector<vector<double> > x2d;
  x2d.resize(10000, vector<double>(3, 1.));
  vector<vector<vector<double> > > x3d;
  x3d.resize(10000, vector<vector<double> >(1, vector<double>(3, 1.)));
  double t = 0;
//  for (int n = 0; n < 3; ++n) {
//    for (int i = 0; i < int(x.size()); ++i) {
//      for (int j = 0; j < int(x.size()); ++j) {
//        //t += x[i]*x[j];
//        t += x2d[i][0]*x2d[j][0];
//      }
//    }
//    std::cout << n << " " << t << std::endl;
//  }
  for (int n = 0; n < 3; ++n) {
    for (int i = 0; i < int(x2d.size()); ++i) {
      for (int k = 0; k < 3; ++k) {
      for (int j = 0; j < int(x2d.size()); ++j) {
      for (int g = 0; g < 3; ++g) {
        //t += x2d[i][k]*x2d[j][g];
        t += x3d[i][0][k]*x3d[j][0][g];
      }}}
    }
    std::cout << n << " " << t << std::endl;
  }
}

/**
 * given range of n, exponent, the number of windows, and overlap of windows, return vector with min and max
 */
vector<vector<int> > nWindow(const int nMolMin,   //!< minimum n
  const int nMolMax,   //!< maximum n
  const double nExp,   //!< exponent of distribution
  const int nWindow,   //!< number of windows
  const int nOverlap   //!< overlap of n between windows
  ) {
  vector<vector<double> > nwind;
  long long tmp = pow(nMolMax, nExp) - pow(nMolMin, nExp);
  const long double dn = tmp/double(nWindow);
  for (int i = 0; i < nWindow; ++i) {
    vector<double> nwintmp;
    if (i == 0) {
      nwintmp.push_back(nMolMin);
    } else {
      nwintmp.push_back(nwind[i-1][1] - nOverlap);
    }
    nwintmp.push_back(pow(pow(nwintmp[0], nExp) + dn, 1./nExp)+nOverlap);
    nwind.push_back(nwintmp);
  }
  vector<vector<int> > nwin;
  for (int i = 0; i < nWindow; ++i) {
    vector<int> nwintmp;
    nwintmp.push_back(nwind[i][0]);
    nwintmp.push_back(nwind[i][1]);
    nwin.push_back(nwintmp);
  }
  nwin.front().front() = nMolMin;
  nwin.back().back() = nMolMax;
  //for (int i = 0; i < nWindow; ++i) cout << nwin[i][0] << " " << nwin[i][1] << " " << nwin[i][1] - nwin[i][0] <<  endl;
  return nwin;
}

/**
 * this function is used to determine range of windows in parallelization
 *  the order parameter of the windows is in the range [mMin, mMax]
 *  the number of windows is nWindow
 *  the bin width of the order parameter is dm
 *
 *  the windows vary in sized based on the grow parameter.
 *  if grow == 0, the windows are roughtly the same size
 *  if grow <  0, the windows decrease in size
 *  if grow >  0, the windows increase in size
 *
 *  For example, if grow > 0, n=3, the windows may look like this
 *
 *  |---------------------| total window
 *  |--|------|-----------| three windows
 *
 */
vector<vector<double> > nWindowGrowth(const double mMin,   //!< minimum m
  const double mMax,   //!< maximum m
  const double grow,   //!< percentage to grow each window
  const int nWindow,   //!< number of windows
  const double dm,     //!< width of m
  const int overlap    //!< overlap+1 overlapping windows
  ) {
  vector<vector<double> > win(nWindow, vector<double>(2));
  if (grow != 0) myOutF("error", std::ostringstream().flush() << "nWindowGrow isn't correctly implemented for grow(" << grow << ") != 0");
  //cout << "mMin " << mMin << endl;

  // first, determine the size of the first window, d0
  //  this value may need to be shifted afterward to fit the bins
  int series = 0;
  for (int i = 1; i < nWindow-1; ++i) series += (nWindow-i);
  //cout << "series " << series << endl;
  double g = grow*(mMax - mMin);
  double d0 = ((mMax - mMin) - g*double(series))/double(nWindow);
  d0 = d0 - fmod(d0, dm);
  //cout << "d0 " << d0 << endl;

  // now that d0 is obtained, the actual grow is computed
  double gactual = (mMax-mMin-double(nWindow)*d0)/double(series);
  gactual = gactual - fmod(gactual, dm);
  if (grow == 0) gactual = 0;
  //cout << "gactual " << gactual << endl;

  for (int i = 0; i < nWindow; ++i) {
    if (i == 0) {
      win[i][0] = mMin;
    } else {
      win[i][0] = win[i-1][1] - (overlap+1)*dm;
    }

    if (i == nWindow - 1) {
      win[i][1] = mMax;
    } else {
      int ser = 0;
      for (int ii = i; ii > 0; --ii) ser += ii;
      win[i][1] = (mMin+(int(overlap/2)+1)*dm) + (i+1)*d0 + gactual*ser;
    }
  }
  //for (int i = 0; i < nWindow; ++i) cout << win[i][0] << " " << win[i][1] << " " << win[i][1] - win[i][0] <<  endl;
  return win;
}
vector<vector<double> > nWindowGrowth(const double mMin, const double mMax, const double grow, const int nWindow, const double dm) {
  return nWindowGrowth(mMin, mMax, grow, nWindow, dm, 0);
}

///**
// * function to return quaternions, given current positions and reference positions
// *  use axis, angle definition of quaternion
// *  qx = ax * sin(angle/2)
// *  qy = ay * sin(angle/2)
// *  qz = az * sin(angle/2)
// *  qw = cos(angle/2)
// *  ax*ax + ay*ay + az*az = 1
// */
//vector<double> pos2quat(const vector<vector<double> > x,   //!< current position
//  const vector<vector<double> > xref    //!< reference position
//  ) {
//
//  // axis is defined as the vector connecting the first and second atoms, normalized to one
//  vector<double> a;
//  double aNorm;
//  for (int dim = 0; dim < int(x[0].size()); ++dim) {
//    a.push_back(x[0][dim] - xref[0][dim]);
//    aNorm += pow(a.back(), 2);
//  }
//  aNorm = sqrt(aNorm);
//  for (int i = 0; i < int(a.size()); ++i) a[i] /= aNorm;
//
//  // angle is defined as the angle
//
//  vector<double> q;
//  return q;
//}
//
///**
// * function to return positions, given current quaternions and reference positions
// */
//vector<vector<double> > quat2pos(const vector<double> q,   //!< current quaternion
//  const vector<vector<double> > xref    //!< reference position
//  ) {
//  vector<vector<double> > x;
//  return x;
//}

/**
 * given a file name in char*, return string which trims all characters proceeding and up to the last specialchr
 */
std::string myTrim(const char* specialchr,  //!< character with which to search for to decide where to trim
  const char* fileName //!< filename to trim
  ) {
  std::string fs(fileName);

  // find position of last spchr
  std::string spchr(specialchr);
  std::size_t found = fs.find(spchr);
  std::size_t foundprev = 0;
  while (found != std::string::npos) {
    foundprev = found;
    found = fs.find(spchr, found + 1);
  }

  // erase all characters up to spchr
  if (int(foundprev) != 0) {
    fs.erase(fs.begin(), fs.begin() + foundprev + spchr.size());
  }

  return fs;
}
std::string myTrim(const char* specialchr, string fileName) {
  return myTrim(specialchr, fileName.c_str());
}

/**
 * funciton to output error messages
 */
void myOut(const char* type,  //!< type of message: verbose, warning or error
  std::ostream& message,  //!< error message
  const std::string className,   //!< name of class experiencing error
  const int verbose
) {
  std::string mstr = dynamic_cast<std::ostringstream&>(message).str();
  std::string typestr(type);
  if ( (verbose == 1) || (typestr.compare("verbose") != 0) ) {
    //const int nproc = omp_get_thread_num();
    int nproc = 0;
    #ifdef MPI_H_
      int initialized;
      MPI_Initialized(&initialized);
      if (initialized) MPI_Comm_rank(MPI_COMM_WORLD, &nproc);
    #endif  // MPI_H_
    #ifdef _OPENMP
      nproc = omp_get_thread_num();
      #pragma omp critical
      {
    #endif  // _OPENMP
    cout << "# " << type << " on proc " << nproc << " of " << className << ": " << mstr << endl;
    #ifdef _OPENMP
      }
    #endif  // _OPENMP
    message.clear();
    if (typestr.compare("error") == 0) exit(0);
  }
}

/**
 * function to print error messages in this particular class
 */
void myOutF(const char* messageType, std::ostream& message) {
  const std::string s("Functions");
  myOut(messageType, message, s, 0);
}

/**
  * function to skip all lines beginning with character in file
 */
void skipCharsInFile(const char comment, std::ifstream &file) {
  std::string line;
  getline(file, line);
  while (line[0] == comment) {
    getline(file, line);
  }
}

/// function to skip all lines until reaching a certain line
void readUntil(const char* searchString, //!< stop when this line is reached
  std::ifstream &file   // file to read
  ) {
  std::string line;
  getline(file, line);
  const int nMax = 1e7;
  int i = 0;
  while (line.compare(searchString) != 0) {
    getline(file, line);
    ++i;
    if (i > nMax) myOutF("error", std::ostringstream().flush() << "readUntil reached nMax attempts(" << nMax << ") while looking for string(" << searchString << ") in file");
    //cout << "readUntil " << i << "/" << nMax << ": " << line << endl;
  }
}

/**
 * return number reported after first appearance of string in file
 */
string fstos(const char* searchString, //!< stop when this line is reached
  const char* fileName // file to read
  ) {
  std::ifstream file(fileName);
  string line;
  string searchStringStr(searchString);
  getline(file, line);
  while ( (line.find(searchString) == std::string::npos) && (file) ) {
    getline(file, line);
  }
  //if (!file) myOutF("error", std::ostringstream().flush() << "can't find searchString(" << searchString << ") in file(" << fileName << ")");
  if (file) {
    line.erase(line.begin(), line.begin() + line.find(searchString) + searchStringStr.size() + 1);
    return line;
  } else {
    string strtmp("");
    return strtmp;
  }
}

/**
 * return number reported after first appearance of string in file
 */
double fstod(const char* searchString, //!< stop when this line is reached
  const char* fileName // file to read
  ) {
  const string str = fstos(searchString, fileName);
  if (str.empty()) {
    myOutF("warning", std::ostringstream().flush() << "can't find searchString(" << searchString << ") in file(" << fileName << ") when converting to double");
    return 0;
  } else {
    return stod(str);
  }
}

/**
 * return number reported after first appearance of string in file
 */
int fstoi(const char* searchString, //!< stop when this line is reached
  const char* fileName // file to read
  ) {
  const string str = fstos(searchString, fileName);
  if (str.empty()) {
    myOutF("warning", std::ostringstream().flush() << "can't find searchString(" << searchString << ") in file(" << fileName << ") when converting to int");
    return 0;
  } else {
    return stoi(str);
  }
}

/**
 * return number reported after first appearance of string in file
 */
unsigned long long fstoull(const char* searchString, //!< stop when this line is reached
  const char* fileName // file to read
  ) {
  const string str = fstos(searchString, fileName);
  if (str.empty()) {
    myOutF("warning", std::ostringstream().flush() << "can't find searchString(" << searchString << ") in file(" << fileName << ") when converting to int");
    return 0;
  } else {
    return stoull(str);
  }
}

/**
 * return number reported after first appearance of string in file
 */
long long fstoll(const char* searchString, //!< stop when this line is reached
  const char* fileName // file to read
  ) {
  const string str = fstos(searchString, fileName);
  if (str.empty()) {
    myOutF("warning", std::ostringstream().flush() << "can't find searchString(" << searchString << ") in file(" << fileName << ") when converting to int");
    return 0;
  } else {
    return stoll(str);
  }
}

/**
 * compute eigenvalues and eigenvectors of a 3x3 real symmetric matrix
 *  based on Jacobi rotations
 *  adapted from LAMMPS, which was adapted from Numerical Recipes jacobi() function
 */
int jacobi(vector<vector<double> > matrix,
  vector<double> &evalues,              //!< eigen values
  vector<vector<double> > &evectors     //!< eigen vectors
  ) {

  int i,j,k;
  double tresh,theta,tau,t,sm,s,h,g,c,b[3],z[3];
  const double MAXJACOBI = 50;

  evalues.resize(matrix.size());
  evectors.resize(matrix.size(), vector<double>(matrix.size()));
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) evectors[i][j] = 0.0;
    evectors[i][i] = 1.0;
  }
  for (i = 0; i < 3; i++) {
    b[i] = evalues[i] = matrix[i][i];
    z[i] = 0.0;
  }

  for (int iter = 1; iter <= MAXJACOBI; iter++) {
    sm = 0.0;
    for (i = 0; i < 2; i++) {
      for (j = i+1; j < 3; j++) {
        sm += fabs(matrix[i][j]);
        if (sm == 0.0) return 0;
      }
    }

    if (iter < 4) {
      tresh = 0.2*sm/(3*3);
    } else {
      tresh = 0.0;
    }

    for (i = 0; i < 2; i++) {
      for (j = i+1; j < 3; j++) {
	      g = 100.0*fabs(matrix[i][j]);
	      if (iter > 4 && fabs(evalues[i])+g == fabs(evalues[i])
	          && fabs(evalues[j])+g == fabs(evalues[j])) {
	        matrix[i][j] = 0.0;
	      } else if (fabs(matrix[i][j]) > tresh) {
	        h = evalues[j]-evalues[i];
	        if (fabs(h)+g == fabs(h)) {
            t = (matrix[i][j])/h;
	        } else {
	          theta = 0.5*h/(matrix[i][j]);
	          t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	          if (theta < 0.0) t = -t;
	        }
          c = 1.0/sqrt(1.0+t*t);
          s = t*c;
          tau = s/(1.0+c);
          h = t*matrix[i][j];
          z[i] -= h;
          z[j] += h;
          evalues[i] -= h;
          evalues[j] += h;
          matrix[i][j] = 0.0;
          for (k = 0; k < i; k++) rotateJacobi(matrix,k,i,k,j,s,tau);
          for (k = i+1; k < j; k++) rotateJacobi(matrix,i,k,k,j,s,tau);
          for (k = j+1; k < 3; k++) rotateJacobi(matrix,i,k,j,k,s,tau);
          for (k = 0; k < 3; k++) rotateJacobi(evectors,k,i,k,j,s,tau);
        }
      }
    }

    for (i = 0; i < 3; i++) {
      evalues[i] = b[i] += z[i];
      z[i] = 0.0;
    }
  }
  myOutF("error", std::ostringstream().flush() << "jacobi failed");
  return 1;
}

/**
 * preform a single Jacobi rotation
 *  adapted from LAMMPS, via NR
 */
void rotateJacobi(vector<vector<double> > &matrix, const int i, const int j, const int k, const int l, const double s, const double tau) {
  const double g = matrix[i][j];
  const double h = matrix[k][l];
  matrix[i][j] = g-s*(h+g*tau);
  matrix[k][l] = h+s*(g-h*tau);
}

/// returns a unit vector orthogonal to given vector
vector<double> orthogonalVec(const vector<double> x) {

  if (int(x.size()) != 3) myOutF("error", std::ostringstream().flush() << "orthogonalVec requires dimensions(" << x.size() << ") == 3");

  // find the component with the least absolute value
  int minIndex = 0;
  if (fabs(x[1]) < fabs(x[minIndex])) minIndex = 1;
  if (fabs(x[2]) < fabs(x[minIndex])) minIndex = 2;

  // define orthogonal vector by setting the minIndex to 0 and swapping the other two indices, setting one of the swapped negative
  vector<double>r(3);
  if (minIndex == 0) {
    r[0] = 0.;
    r[1] = -x[2];
    r[2] = x[1];
  } else if (minIndex == 1) {
    r[0] = -x[2];
    r[1] = 0.;
    r[2] = x[0];
  } else {
    r[0] = -x[1];
    r[1] = x[0];
    r[2] = 0.;
  }

  // normalize
  const double d = sqrt(myVecDotProd(r, r));
  for (int dim = 0; dim < 3; ++dim) {
    r[dim] /= d;
  }

  vector<double> x2 = x;
  if (fabs(myVecDotProd(r, x2)) > doubleTolerance) myOutF("error", std::ostringstream().flush() << "not orthogonal");

  return r;
}

/// returns rotation matrix for arbitrary rotation axis, u, by an angle theta
vector<vector<double> > rotMatAxisAngle(vector<double> u, const double theta) {
  if (int(u.size()) != 3) myOutF("error", std::ostringstream().flush() << "rotMatAxisAngle requires dimensions(" << u.size() << ") == 3");
  if (fabs(myVecDotProd(u, u) - 1) > doubleTolerance) myOutF("error", std::ostringstream().flush() << "axis u(" << u[0] << "," << u[1] << "," << u[2] << ") in rotMatAxisAngle must be normalized");
  vector<vector<double> > r(3, vector<double>(3));
  const double c=cos(theta), s=sin(theta);
  r[0][0] = u[0]*u[0]*(1-c)+c;
  r[0][1] = u[0]*u[1]*(1-c)-u[2]*s;
  r[0][2] = u[0]*u[2]*(1-c)+u[1]*s;
  r[1][0] = u[0]*u[1]*(1-c)+u[2]*s;
  r[1][1] = u[1]*u[1]*(1-c)+c;
  r[1][2] = u[1]*u[2]*(1-c)-u[0]*s;
  r[2][0] = u[0]*u[2]*(1-c)-u[1]*s;
  r[2][1] = u[1]*u[2]*(1-c)+u[0]*s;
  r[2][2] = u[2]*u[2]*(1-c)+c;
//  for (int i = 0; i < 3; ++i) {
//    for (int j = 0; j < 3; ++j) {
//      cout << r[i][j] << " ";
//    }
//    cout << endl;
//  }
  return r;
}

/// returns vector rotated by angle theta about axis u
vector<double> rotateVecByAxisAngle(vector<double> x, vector<double> axis, const double theta) {
  vector<vector<double> > r = rotMatAxisAngle(axis, theta);
  return myMatVecMul(r, x);
}

/// draw a superquadric
//  |x/A|^r + |y/B|^s + |z/C|^t <= 1
vector<vector<double> > superQuadric3D(const double A, const double m, const int n, const double umid, const double vmid, const double du, const double dv) {
  vector<vector<double> > xAll;
  vector<double> x(3);
  const double ddvv = 2*dv/double(n-1), dduu = 2*du/double(n-1);
  for (double v = vmid-dv; v <= vmid+dv; v+=ddvv) {
    for (double u = umid-du; u <= umid+du; u+=dduu) {
//        cout << "v " << v << " u " << u << " cv " << cos(v) << " sv " << sin(v) << " cu " << cos(u) << " su " << sin(u) << endl;
      x[0] = A*cosfSuperQuad(v, 2./m)*cosfSuperQuad(u, 2./m);
//        cout << "cosfsq " << cosfSuperQuad(v, 2./m)<< endl;
      x[1] = A*cosfSuperQuad(v, 2./m)*sinfSuperQuad(u, 2./m);
      x[2] = A*sinfSuperQuad(v, 2./m);
//        cout << "A " << A[0] << " " << A[1] << " " << A[2] << endl;
//        cout << "eps " << m << " " << m << " " << m << endl;
//        cout << x[0] << " " << x[1] << " " << x[2] << endl;
      xAll.push_back(x);
    }
  }
  return xAll;
}
vector<vector<double> > superQuadric3D(const double A, const double m, const int n) {
  return superQuadric3D(A, m, n, 0., 0., PI, PI/2.);
}
double cosfSuperQuad(const double w, const double m) {
  return mySign(pow(fabs(cos(w)), m), cos(w));
}
double sinfSuperQuad(const double w, const double m) {
  return mySign(pow(fabs(sin(w)), m), sin(w));
}

/**
 * minimum distance between two objects described by grid points
 */
void minDistGrid(const vector<vector<double> > &x1, const vector<vector<double> > &x2, double &xmin, int &x1min, int &x2min) {
  if (int(x1[0].size()) != 3) myOutF("error", std::ostringstream().flush() << "mindistGrid assume 3d");
  double xSqMin = 1e100;
  for (int i1 = 0; i1 < int(x1.size()); ++i1) {
    const double xi = x1[i1][0],
                 yi = x1[i1][1],
                 zi = x1[i1][2];
    for (int i2 = 0; i2 < int(x1.size()); ++i2) {
      const double dx = xi - x2[i2][0],
                   dy = yi - x2[i2][1],
                   dz = zi - x2[i2][1];
      const double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < xSqMin) {
        x1min = i1;
        x2min = i2;
        xSqMin = r2;
      }
    }
  }
  xmin = sqrt(xSqMin);
}

// convert cartesian coordinates vector to spherical coordinates
vector<double> cartesian2spherical(vector<double> rCart) {

  if (int(rCart.size()) != 3) myOutF("error", std::ostringstream().flush() << "cartesian2spherecal is only implemented for dimensionality (" << rCart.size() << ") of 3");
  vector<double> rSphere(rCart.size());
  rSphere[0] = sqrt(myVecDotProd(rCart, rCart));
  rSphere[1] = acos(rCart[2]/rSphere[0]);
  rSphere[2] = atan2(rCart[1], rCart[0]);
  vector<double> rCartCheck(rCart.size());
  rCartCheck[0] = rSphere[0]*sin(rSphere[1])*cos(rSphere[2]);
  rCartCheck[1] = rSphere[0]*sin(rSphere[1])*sin(rSphere[2]);
  rCartCheck[2] = rSphere[0]*cos(rSphere[1]);
  const double check = fabs(myVecDotProd(rCart, rCart) - myVecDotProd(rCartCheck, rCartCheck));
  if (check > doubleTolerance) myOutF("error", std::ostringstream().flush() << "cartesian2spherical check failure " << check);
  return rSphere;
}

// returns the spherical harmonics (l=6) of a cartesian vector r
vector<std::complex<double> > cart2sphereHarm6(vector<double> rCart) {
  vector<double> rs = cartesian2spherical(rCart);
  vector<double> sphHr(13);
  vector<std::complex<double> > sphH(13);
  const double theta = rs[1], phi = rs[2], st = sin(theta), ct = cos(theta);

  sphHr[0] = 1./64.*sqrt(3003./PI)*pow(st, 6);
  sphHr[1] = -3./32.*sqrt(1001./PI)*pow(st, 5)*ct;
  sphHr[2] = 3./32.*sqrt(91./2./PI)*pow(st, 4)*(11*pow(ct,2) - 1.);
	sphHr[3] = -1./32.*sqrt(1365./PI)*pow(st,3)*(11.*pow(ct,3) - 3.*ct);
	sphHr[4] = 1./64.*sqrt(1365./PI)*pow(st,2)*(33.*pow(ct,4) - 18.*pow(ct,2) + 1.);
	sphHr[5] =-1./16.*sqrt(273./2./PI)*st;
	sphHr[5] *= 33.*pow(ct,5) - 30.*pow(ct,3) + 5.*ct;
	sphHr[6] = 1./32.*sqrt(13./PI)*(231.*pow(ct,6) - 315.*pow(ct,4) + 105.*pow(ct,2) - 5.);
	sphHr[7] = sphHr[5]*-1.;
	sphHr[8] = sphHr[4];
	sphHr[9] = sphHr[3]*-1;
	sphHr[10] = sphHr[2];
  sphHr[11] = sphHr[1]*-1;
	sphHr[12] = sphHr[0];

	for (int i = 0; i < 13; ++i) {
    const int l = 6-i;
    sphH[i] = sphHr[i] * std::complex<double>(cos(double(l)*phi), sin(double(l)*phi));
  }
  return sphH;
}

// returns sum of vector of complex numbers multiplied by its conjugate
double complexVec2norm(vector<std::complex<double> > compVec) {
  double norm = 0.;
  for (int i = 0; i < int(compVec.size()); ++i) {
    norm += (compVec[i] * std::conj(compVec[i])).real();
  }
  return norm;
}

/**
 * round double to nearest integer
 */
int myRound(double x) {
  return floor(x + 0.5);
}

/**
 * generate a rotation matrix from euler angles
 *
 * http://mathworld.wolfram.com/EulerAngles.html
 *
 * Euler angles (phi, theta, psi) are defined by rotation phi about z-axis, rotation theta about new x-xis, then rotation psi about new z-axis
 *  phi [-pi,pi], theta [0, pi], psi [-pi, pi]
 */
vector<vector<double> > Euler2RotMat(const vector<double> euler) {
  if (euler.size() != 3) myOutF("error", std::ostringstream().flush() << "assumes 3D, but size of euler is " << euler.size());
  vector<vector<double> > R(3, vector<double>(3));
  const double sphi = sin(euler[0]),
               cphi = cos(euler[0]),
               stheta = sin(euler[1]),
               ctheta = cos(euler[1]),
               spsi = sin(euler[2]),
               cpsi = cos(euler[2]);
  R[0][0] = cpsi*cphi - ctheta*sphi*spsi; //a11
  R[0][1] = cpsi*sphi + ctheta*cphi*spsi; //a12
  R[0][2] = spsi*stheta;                  //a13
  R[1][0] = -spsi*cphi - ctheta*sphi*cpsi;//a21
  R[1][1] = -spsi*sphi + ctheta*cphi*cpsi;//a22
  R[1][2] = cpsi*stheta;                  //a23
  R[2][0] = stheta*sphi;                  //a31
  R[2][1] = -stheta*cphi;                 //a32
  R[2][2] = ctheta;                       //a33
  return R;
}

/**
 * obtain euler angles from a rotation matrix (same convention as Euler2RotMat)
 */
vector<vector<double> > RotMat2Euler(const vector<vector<double> > rotMat) {
  if (rotMat.size() != 3) myOutF("error", std::ostringstream().flush() << "assumes 3D, but size of rotMat is " << rotMat.size());
  if (rotMat[0].size() != 3) myOutF("error", std::ostringstream().flush() << "assumes 3D, but size of rotMat[0] is " << rotMat[0].size());
  double theta1, phi1, psi1;
  const double ctheta1 = rotMat[2][2];
  if ( fabs(ctheta1-1.) < 1e-12) {
    theta1 = acos(1.);
  } else {
    theta1 = acos(ctheta1);
  }
  const double stheta = sin(theta1);
  //cout << "theta1 " << theta1 << " stheta " << stheta << endl;
  if (fabs(stheta) < 1e-8) {
    phi1 = 0;
    psi1 = atan2(rotMat[0][1], rotMat[0][0]);
  } else {
    phi1 = atan2(rotMat[2][0]/stheta,-rotMat[2][1]/stheta);
    psi1 = atan2(rotMat[0][2]/stheta, rotMat[1][2]/stheta);
  }
  //cout << "phi1 " << phi1 << " psi1 " << psi1 << endl;
  vector<vector<double> > euler(2, vector<double>(3));
  euler[0][0] = phi1;
  euler[0][1] = theta1;
  euler[0][2] = psi1;
  return euler;
}

/**
 * determinant of a 3x3 matrix
 */
double det3by3(const vector<vector<double> > mat) {
  if (mat.size() != 3) myOutF("error", std::ostringstream().flush() << "assumes 3D, but size of mat is " << mat.size());
  if (mat[0].size() != 3) myOutF("error", std::ostringstream().flush() << "assumes 3D, but size of mat[0] is " << mat[0].size());
  const double a = mat[0][0], b = mat[0][1], c = mat[0][2], d=mat[1][0], e=mat[1][1], f=mat[1][2], g=mat[2][0], h=mat[2][1], i=mat[2][2];
  return a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;
}

/**
 * determinant of a 2x2 matrix
 */
double det2by2(const vector<vector<double> > mat) {
  if (mat.size() != 2) myOutF("error", std::ostringstream().flush() << "assumes 2D, but size of mat is " << mat.size());
  if (mat[0].size() != 2) myOutF("error", std::ostringstream().flush() << "assumes 2D, but size of mat[0] is " << mat[0].size());
  const double a = mat[0][0], b = mat[0][1], c = mat[1][0], d=mat[0][1];
  return a*d-b*c;
}

/**
 * cofactor of a 3x3 matrix
 */
vector<vector<double> > cofactor3by3(const vector<vector<double> > mat) {
  if (mat.size() != 3) myOutF("error", std::ostringstream().flush() << "assumes 3D, but size of mat is " << mat.size());
  if (mat[0].size() != 3) myOutF("error", std::ostringstream().flush() << "assumes 3D, but size of mat[0] is " << mat[0].size());
  vector<vector<double> > cof(3, vector<double>(3));
  cof[0][0] =  mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1];
  cof[0][1] = -mat[1][0]*mat[2][2] + mat[1][2]*mat[2][0];
  cof[0][2] =  mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0];

  cof[1][0] = -mat[0][1]*mat[2][2] + mat[0][2]*mat[2][1];
  cof[1][1] =  mat[0][0]*mat[2][2] - mat[0][2]*mat[2][0];
  cof[1][2] = -mat[0][0]*mat[2][1] + mat[0][1]*mat[2][0];

  cof[2][0] =  mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1];
  cof[2][1] = -mat[0][0]*mat[1][2] + mat[0][2]*mat[1][0];
  cof[2][2] =  mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
  return cof;
}

/**
 * transpose
 */
vector<vector<double> > transpose(const vector<vector<double> > mat) {
  vector<vector<double> > transpose(int(mat[0].size()), vector<double>(int(mat.size())));
  for (int idim = 0; idim < int(mat.size()); ++idim) {
    for (int jdim = 0; jdim < int(mat[idim].size()); ++jdim) {
      transpose[jdim][idim] = mat[idim][jdim];
    }
  }
  return transpose;
}

/*
 * inverse
 */
vector<vector<double> > inv3by3(const vector<vector<double> > mat) {
  vector<vector<double> > inv = transpose(cofactor3by3(mat));
  const double det = det3by3(mat);
  for (int idim = 0; idim < int(inv.size()); ++idim) {
    for (int jdim = 0; jdim < int(inv[idim].size()); ++jdim) {
      inv[idim][jdim] /= det;
    }
  }
  return inv;
}

/// convert quaternion vector to euler angles (uses quat2rot -> RotMat2euler)
vector<double> quat2euler(vector<double> quat) {
  return RotMat2Euler( (quat2rot(quat)) )[0];
  //return RotMat2Euler( transpose(quat2rot(quat)) )[0];
}

/// convert rotation matrix to quaternions
vector<double> rot2quat(const vector<vector<double> > &rot) {
  vector<double> q(4);

  // begin by attempting to invert equations for quat2rot
  //  manipulate rotMat to find two equations involing q0 and q1
  //  r10+r01 = 4q0q1 = A1
  //  r00-r11 = 2q0^2 - 2q1^2 = A2
  //  A2q1^2 = 2A1^2 - 2q1^4
  //  for x = q1^2,
  //    (-2)x^2 + (-A2)x + (2A1^2) = 0
  //    assuming q1 != 0, so check discriminant of quadratice equation

  double x1, x2;
  const double discrim = quadraticEqReal(-2., -(rot[0][0]-rot[1][1]), 2*pow(rot[1][0]+rot[0][1], 2), x1, x2);
  if (discrim <= 0) {
    cout << "oh no!" << endl;
    exit(0);
  }

  if (x1 > 0) {
    q[1] = sqrt(x1);
  } if (x2 > 0) {
    q[1] = sqrt(x2);
  } else {
    cout << "not again" << endl;
    exit(0);
  }

  // asuming q1 != 0...
  q[0] = (rot[1][0]+rot[0][1])/4./q[1];
  q[2] = (rot[2][0]+rot[0][2])/4./q[0];
  q[3] = sqrt( (rot[0][0]+rot[1][1]+2*pow(q[2],2))/2);
  return q;
}

/**
 * periodic boundary conditions for euler angles
 *
 * Euler angles (phi, theta, psi) are defined by rotation phi about z-axis, rotation theta about new x-xis, then rotation psi about new z-axis
 *  phi [-pi,pi], theta [0, pi], psi [-pi, pi]
 */
void eulerPBC(vector<double> &euler) {
  while (euler[0] >  PI) euler[0] -= 2*PI;
  while (euler[0] < -PI) euler[0] += 2*PI;
  while (euler[1] >  PI) euler[1] -= PI;
  while (euler[1] <  0)  euler[1] += PI;
  while (euler[2] >  PI) euler[2] -= 2*PI;
  while (euler[2] < -PI) euler[2] += 2*PI;
}

/**
 * my acos includes some double precision checks
 */
double myACOS(double x) {
  if (fabs(x-1.) < doubleTolerance) {
    return acos(1.);
  } else if (fabs(x+1.) < doubleTolerance) {
    return acos(-1.);
  } else {
    return acos(x);
  }
}

