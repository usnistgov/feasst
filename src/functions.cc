/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "functions.h"
#include <fstream>
#include <algorithm>
#include <math.h>
#include <complex>
#include <sys/stat.h>
#include <limits>
#include <signal.h>

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

double sign(const double a,	const double b) {
  return (b >= 0.0 ? fabs(a) : -fabs(a));
}

/// Function to fill 2D vector with input variable.
void fill(const double input, vector<vector<double> > &x) {
  for (vector<vector<double> >::iterator iter = x.begin();
       iter != x.end(); ++iter) {
    std::fill((*iter).begin(), (*iter).end(), input);
  }
}

void fill(const int input, vector<vector<int> > &x) {
  for (vector<vector<int> >::iterator iter = x.begin();
       iter != x.end(); ++iter) {
    std::fill((*iter).begin(), (*iter).end(), input);
  }
}

int numLines(const std::string fileName) {
  std::ifstream file(fileName.c_str());
  ASSERT(file.good(), "cannot open file " << fileName);

	int n = 0;
	std::string line;
	while(std::getline(file, line)) ++n;
	return n;
}

vector<double> matVecMul(const vector<vector<double> > &a,
	const vector<double> &x) {
  vector<double> b(int(a.size()));
	ASSERT(static_cast<int>(a[0].size()) == static_cast<int>(x.size()),
    "matVecMul requires columns of matrix a(" << a[0].size()
    << ") to equal rows of vector x(" << x.size());
	for (int j = 0; j < int(a.size()); ++j) {
		for (int i = 0; i < int(x.size()); ++i) {
			b[j] += a[j][i] * x[i];
		}
	}

  return b;
}

double vecDotProd(const vector<double> &a, const vector<double> &b) {
  ASSERT(a.size() == b.size(), "vecDotProd requires vectors of equal size, "
    << "however, a(" << a.size() << ") and b(" << b.size() << ")");
	double c = 0;
	for (int i = 0; i < static_cast<int>(a.size()); ++i) {
		c += a[i]*b[i];
	}
	return c;
}

void normalizeVec(vector<double> *x) {
  const double len = sqrt(vecDotProd(*x, *x));
  for (int dim = 0; dim < int((*x).size()); ++dim) {
    (*x)[dim] /= len;
  }
}

vector<double> crossProd(const vector<double> &a, const vector<double> &b) {
  ASSERT( (a.size() == b.size()) && (a.size() == 3),
    "crossProd requires vectors of 3 dimensions, however, a(" << a.size()
    << ") and b(" << b.size() << ")");
  vector<double> c(a.size());
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
  return c;
}

void ranInitByDate() {
	const int t = time(NULL);
  srand ( t );
  NOTE("time(seed): " << t);
}

void ranInitForRepro(const int seed) {
	srand ( seed );
  NOTE("Initializing random number generator for reproduction with seed("
    << seed << ")");
}

vector<vector<double> > quat2rot(vector<double> q) {
  vector<vector<double> > r(3, vector<double>(3));
  ASSERT(q.size() == 4, "1quaterion is a 4d vector in 3d");
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

vector<vector<double> > theta2rot(double theta) {
  vector<vector<double> > r(2, vector<double>(2));
  r[0][0] = cos(theta);
  r[0][1] = -sin(theta);
  r[1][0] = -r[0][1];
  r[1][1] = r[0][0];
  return r;
}

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

vector<vector<double> > matMul(const vector<vector<double> > &a,
  const vector<vector<double> > &b) {
  vector<vector<double> > c(int(a.size()), vector<double>(int(b[0].size())));
  ASSERT(a[0].size() == b.size(), "matMul requires columns of matrix a("
    << a[0].size() << ") to equal rows of vector x(" << b.size());
  ASSERT(a.size() == c.size(), "matMul requires rows of matrix a("
    << a.size() << ") to equal rows of vector b(" << c.size());
	for (int i = 0; i < int(a.size()); ++i) {
		for (int j = 0; j < int(b[0].size()); ++j) {
		  for (int k = 0; k < int(b.size()); ++k) {
		    c[i][j] += a[i][k]*b[k][j];	
		  }
    }
	}
  return c;
}

double volShell(const double rabove, const double rbelow, const int dim) {
  ASSERT(dim == 3, "only implemented for 3 dimensions (dim=" << dim << ")");
  return 4./3.*PI*(rabove*rabove*rabove - rbelow*rbelow*rbelow);
}

bool fileExists(const char* fileName) {
  struct stat buf;
  if (stat(fileName, &buf) != -1) return true;
  return false;
}
bool fileExists(std::ifstream& file) {
  if (file.peek() == std::ifstream::traits_type::eof()) return false;
  return true;
}

void fileBackUp(const char* fileName) {
  if (fileExists(fileName)) {
    std::ostringstream f;
    f << fileName << ".bak";
    rename(fileName, f.str().c_str());
  }
}
void fileBackUp(const std::string fileName) {
  fileBackUp(fileName.c_str());
}

void skipCharsInFile(const char comment, std::ifstream &file) {
  std::string line;
  getline(file, line);
  while (line[0] == comment) {
    getline(file, line);
  }
}

int readUntil(const char* searchString, std::ifstream &file,
  const int optional) {
  std::string line;
  getline(file, line);
  const int nMax = 1e7;
  int i = 0;
  while ((line.compare(searchString) != 0) && (i != nMax)) {
    getline(file, line);
    ++i;
  }
  // check if not found
  if (i == nMax) {
    if (optional == 0) {
      ASSERT(0, "readUntil reached nMax attempts(" << nMax
        << ") while looking for string(" << searchString << ") in file");
    }
    return 0;
  }
  return 1;
}

vector<vector<int> > nWindow(const int nMolMin, const int nMolMax,
  const double nExp, const int nWindow, const int nOverlap) {
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

vector<vector<double> > nWindowGrowth(const double mMin, const double mMax,
  const double grow, const int nWindow, const double dm, const int overlap) {
  vector<vector<double> > win(nWindow, vector<double>(2));
  ASSERT(grow == 0, "nWindowGrow isn't correctly implemented for grow("
    << grow << ") != 0");

  // first, determine the size of the first window, d0
  //  this value may need to be shifted afterward to fit the bins
  int series = 0;
  for (int i = 1; i < nWindow-1; ++i) series += (nWindow-i);
  double g = grow*(mMax - mMin);
  double d0 = ((mMax - mMin) - g*double(series))/double(nWindow);
  d0 = d0 - fmod(d0, dm);

  // now that d0 is obtained, the actual grow is computed
  double gactual = (mMax-mMin-double(nWindow)*d0)/double(series);
  gactual = gactual - fmod(gactual, dm);
  if (grow == 0) gactual = 0;

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
  return win;
}

std::string trim(const char* specialchr, const char* fileName, int fromLeft) {
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
    if (fromLeft == 1) {
      fs.erase(fs.begin(), fs.begin() + foundprev + spchr.size());
    } else {
      fs.erase(fs.begin() + foundprev + spchr.size(), fs.end());
    }
  }

  return fs;
}
std::string trim(const char* specialchr, string fileName, int fromLeft) {
  return trim(specialchr, fileName.c_str(), fromLeft);
}

string fstos(const char* searchString, const char* fileName) {
  std::ifstream file(fileName);
  string line;
  string searchStringStr(searchString);
  getline(file, line);
  while ( (line.find(searchString) == std::string::npos) && (file) ) {
    getline(file, line);
  }
  if (file) {
    line.erase(line.begin(), line.begin() + line.find(searchString) + searchStringStr.size() + 1);
    return line;
  } else {
    string strtmp("");
    return strtmp;
  }
}

double fstod(const char* searchString, const char* fileName) {
  const string str = fstos(searchString, fileName);
  if (str.empty()) {
    NOTE("can't find searchString(" << searchString << ") in file("
      << fileName << ") when converting to double");
    return 0;
  } else {
    return stod(str);
  }
}

int fstoi(const char* searchString, const char* fileName) {
  const string str = fstos(searchString, fileName);
  if (str.empty()) {
    NOTE("can't find searchString(" << searchString << ") in file("
      << fileName << ") when converting to double");
    return 0;
  } else {
    return stoi(str);
  }
}

long long fstoll(const char* searchString, const char* fileName) {
  const string str = fstos(searchString, fileName);
  if (str.empty()) {
    NOTE("can't find searchString(" << searchString << ") in file("
      << fileName << ") when converting to double");
    return 0;
  } else {
    return stoll(str);
  }
}

unsigned long long fstoull(const char* searchString, const char* fileName) {
  const string str = fstos(searchString, fileName);
  if (str.empty()) {
    NOTE("can't find searchString(" << searchString << ") in file("
      << fileName << ") when converting to double");
    return 0;
  } else {
    return stoull(str);
  }
}

int jacobi(vector<vector<double> > matrix, vector<double> &evalues,
  vector<vector<double> > &evectors) {
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
  ASSERT(0, "jacobi failed");
  return 1;
}

void rotateJacobi(vector<vector<double> > &matrix, const int i, const int j,
  const int k, const int l, const double s, const double tau) {
  const double g = matrix[i][j];
  const double h = matrix[k][l];
  matrix[i][j] = g-s*(h+g*tau);
  matrix[k][l] = h+s*(g-h*tau);
}

vector<double> orthogonalVec(const vector<double> &x) {
  ASSERT(x.size() == 3, "orthogonalVec requires dimensions("
    << x.size() << ") == 3");

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
  const double d = sqrt(vecDotProd(r, r));
  for (int dim = 0; dim < 3; ++dim) {
    r[dim] /= d;
  }
  ASSERT(fabs(vecDotProd(r, x)) < DTOL, "not orthogonal");
  return r;
}

vector<vector<double> > rotMatAxisAngle(vector<double> u, const double theta) {
  ASSERT(u.size() == 3, "rotMatAxisAngle requires dimensions("
    << u.size() << ") == 3");
  ASSERT(fabs(vecDotProd(u, u) - 1) < DTOL, "axis u(" << u[0] << ","
    << u[1] << "," << u[2] << ") in rotMatAxisAngle must be normalized");
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
  return r;
}

vector<double> rotateVecByAxisAngle(vector<double> x, vector<double> axis,
  const double theta) {
  vector<vector<double> > r = rotMatAxisAngle(axis, theta);
  return matVecMul(r, x);
}

vector<double> cartesian2spherical(vector<double> rCart) {
  ASSERT(rCart.size() == 3, "cartesian2spherecal is only implemented for "
    << "dimensionality (" << rCart.size() << ") of 3");
  vector<double> rSphere(rCart.size());
  rSphere[0] = sqrt(vecDotProd(rCart, rCart));
  rSphere[1] = acos(rCart[2]/rSphere[0]);
  rSphere[2] = atan2(rCart[1], rCart[0]);
  vector<double> rCartCheck(rCart.size());
  rCartCheck[0] = rSphere[0]*sin(rSphere[1])*cos(rSphere[2]);
  rCartCheck[1] = rSphere[0]*sin(rSphere[1])*sin(rSphere[2]);
  rCartCheck[2] = rSphere[0]*cos(rSphere[1]);
  const double check = fabs(vecDotProd(rCart, rCart) - vecDotProd(rCartCheck, rCartCheck));
  ASSERT(check < DTOL, "cartesian2spherical check failure " << check);
  return rSphere;
}

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
    const int l = 6 - i;
    sphH[i] = sphHr[i] * std::complex<double>
      (cos(double(l)*phi), sin(double(l)*phi));
  }
  return sphH;
}

double complexVec2norm(vector<std::complex<double> > compVec) {
  double norm = 0.;
  for (int i = 0; i < int(compVec.size()); ++i) {
    norm += (compVec[i] * std::conj(compVec[i])).real();
  }
  return norm;
}

int feasstRound(double x) { return floor(x + 0.5); }

vector<vector<double> > Euler2RotMat(const vector<double> euler) {
  ASSERT(euler.size() == 3, "assumes 3D, but size of euler is "
    << euler.size());
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

vector<vector<double> > RotMat2Euler(const vector<vector<double> > rotMat) {
  ASSERT(rotMat.size() == 3, "assumes 3D, but size of rotMat is "
    << rotMat.size());
  ASSERT(rotMat[0].size() == 3, "assumes 3D, but size of rotMat[0] is "
    << rotMat[0].size());
  double theta1, phi1, psi1;
  const double ctheta1 = rotMat[2][2];
  if ( fabs(ctheta1-1.) < 1e-12) {
    theta1 = acos(1.);
  } else {
    theta1 = acos(ctheta1);
  }
  const double stheta = sin(theta1);
  if (fabs(stheta) < 1e-8) {
    phi1 = 0;
    psi1 = atan2(rotMat[0][1], rotMat[0][0]);
  } else {
    phi1 = atan2(rotMat[2][0]/stheta,-rotMat[2][1]/stheta);
    psi1 = atan2(rotMat[0][2]/stheta, rotMat[1][2]/stheta);
  }
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
  ASSERT(mat.size() == 3, "assumes 3D, but size of mat is " << mat.size());
  ASSERT(mat[0].size() == 3, "assumes 3D, but size of mat[0] is "
    << mat[0].size());
  const double a = mat[0][0], b = mat[0][1], c = mat[0][2], d=mat[1][0],
    e=mat[1][1], f=mat[1][2], g=mat[2][0], h=mat[2][1], i=mat[2][2];
  return a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h;
}

double det2by2(const vector<vector<double> > mat) {
  ASSERT(mat.size() == 2, "assumes 3D, but size of mat is " << mat.size());
  ASSERT(mat[0].size() == 2, "assumes 3D, but size of mat[0] is "
    << mat[0].size());
  const double a = mat[0][0], b = mat[0][1], c = mat[1][0], d=mat[0][1];
  return a*d-b*c;
}

vector<vector<double> > cofactor3by3(const vector<vector<double> > mat) {
  ASSERT(mat.size() == 3, "assumes 3D, but size of mat is " << mat.size());
  ASSERT(mat[0].size() == 3, "assumes 3D, but size of mat[0] is "
    << mat[0].size());
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

vector<double> quat2euler(vector<double> quat) {
  return RotMat2Euler( (quat2rot(quat)) )[0];
}

bool stringInString(const std::string searchString,
  const std::string stringToSearch) {
  size_t pos;
  pos = stringToSearch.find(searchString);
  if (pos != std::string::npos) return true;
  return false;
}

std::shared_ptr<std::ifstream> make_ifstream(const char* fileName) {
  return std::make_shared<std::ifstream>(fileName);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_
