/**
 *
 * \file
 *
 * \brief external function library
 *
 * Utility functions for multidimensional vectors
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <limits>
#include <memory>
#include <deque>
#include <complex>
#ifdef MPI_H_
  #include <mpi.h>
#endif
#ifdef OMP_H_
  #include <omp.h>
#endif  // OMP_H_
#include "./custom_exception.h"
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::shared_ptr;
using std::make_shared;

# define ASSERT(condition, message) \
if (! (condition)) { \
  std::stringstream err_msg; \
  err_msg << "Assertion `" #condition "` failed in " << __FILE__ \
            << " line " << __LINE__ << ": " << message; \
  customException c(err_msg); \
}
//throw c;

/// function to return the magnitude of the first argument, with the sign of the
//  second argument
double mySign(const double a, const double b);

/// function to add multidimensional vectors
void myAdd(vector<vector<double> > &x, vector<vector<double> > &y);

/// function to fill 2-d vector with input variable
void myFill(const double y, vector<vector<double> > &x);
void myFill(const int y, vector<vector<int> > &x);

/// function to fill 3-d vector with input variable
template<class T>
void myFill(const T y, vector<vector<vector<T> > > &x) {
  for (typename vector<vector<vector<T> > >::iterator iter = x.begin();
       iter != x.end(); ++iter) {
    for (typename vector<vector<T> >::iterator iter2 = (*iter).begin();
         iter2 != (*iter).end(); ++iter2) {
      std::fill((*iter2).begin(), (*iter2).end(), y);
    }
  }
};

/// function to return determinant of 3d symmetric matrix given by 6d vector
double det3DSym(const vector<double> &x);

/// function to return trace of 3d symmetric matrix given by 6d vector
double tr3DSym(const vector<double> &x);

/// function to return square of two 3d symmetric matrix given by 6d vector
void sq3DSym(vector<double> x, vector<double> &y);

/// function to return the real roots of cubic equation given by 4d vector
void cubicSolve(const double b, const double c, const double d,
                double &r1, double &r2, double &r3);
//void cubicSolve(const double a, const double b, const double c, const double d, vector<double> &y);
//void cubicSolve(vector<double> &a, vector<double> &y);

/// function to return the number of lines in a file
int numLines(const string fileName);

/// function to return vector which is product of matrix and vector
void myMatVecMul(vector<vector<double> > &a, vector<double> &x, vector<double> &b);
vector<double> myMatVecMul(vector<vector<double> > &a, vector<double> &x);

/// function to return vector dot product, a.b=c
double myVecDotProd(vector<double> &a, vector<double> &b);

/// function to normalize vector to unit size
void myVecNormalize(vector<double> &x);

/// function to return vector cross product, a*b=c
void myVecCrosProd(vector<double> &a, vector<double> &b, vector<double> &c);

/// minimum image separation vector
void minimumImage(const vector<double> &xi, const vector<double> &xj, const vector<vector<double> > &boxMat, vector<double> &xij);

/// generate spherical grid
void sphereGrid(const int nTheta, const int nPhi, vector<vector<double> > &x);

/// PI = 3.14159265358979323846264338327950 truncated to double precision
const double PI = 4.0*atan(1.0);

/// constants from CODATA
// http://dx.doi.org/10.1103/RevModPhys.84.1527
// http://dx.doi.org/10.1063/1.4724320
const double boltzmannConstant = 1.3806488E-23; //J/K
const double avogadroConstant = 6.02214129E+23; //1/mol
const double elementaryCharge = 1.602176565E-19;  //C
const double permitivityVacuum = 8.854187817E-12;  //C^2/J/m
const double idealGasConstant = boltzmannConstant*avogadroConstant;
const double joulesPercal = 4.184;
const double doubleTolerance = 1e-15;
const double DTOL = doubleTolerance;  // FIX: depreciate doubleTolerance

/// function to return random number (0,1) which requires initialization
void myRanInitByDate();
void myRanInitForRepro();
void myRanInitForRepro(const int seed);
/// function to randomly pick an integer number between bounds min and max, inclusive

/// function to generate rotation matrix from quaternions
void quat2rot(vector<double> &q, vector<vector<double> > &r);
vector<vector<double> > quat2rot(vector<double> q);

/// function to generate 2d rotation matrix from angle
vector<vector<double> > theta2rot(double theta);

/// function to multiply all elements of a one dimensional vector
template<class T>
T myProd(const vector<T> &vec) {
  T prod = 1;
  for (int i = 0; i < int(vec.size()); ++i) {
    prod *= vec[i];
  }
  return prod;
};

/// function to return vector x with atomic positions of a reference SPC/E water molecule
vector<vector<double> > vecSPCE();

/// function to return vector x with atomic positions of a reference one-patch molecule
vector<vector<double> > vecOnePatch();

///// function to return quaternions, given current positions and reference positions
//vector<double> pos2quat(vector<vector<double> > x, vector<vector<double> > xref);

/// function to return product of two matrices
void myMatMul(vector<vector<double> > &a, vector<vector<double> > &b, vector<vector<double> > &c);
vector<vector<double> > myMatMul(vector<vector<double> > &a, vector<vector<double> > &b);

/// function to return complimentary error function, from Numerical Recipes 6.2 pg 221
double erfcc(const double x);

double volShell(const double rabove, const double rbelow, const int dim);   //!< volume of shell

/// function to return average of vector
double myVecAv(const vector<int> &x);
template<class T>
double myVecAv(const vector<T> x) {return std::accumulate(x.begin(), x.end(), 0.) / int(x.size()); };

/// function to open file to write and create back-up
void fileBackUp(const char* fileName);
bool myFileExists(const char* fileName);

/// function to skip all lines beginning with character in file
void skipCharsInFile(const char comment, std::ifstream &file);

/// function to skip all lines until reaching a certain line
void readUntil(const char* searchString, std::ifstream &file);

/// test the speed of arrays
void arraySpeedTest();

/// given range of n, exponent, the number of windows, and overlap of windows, return vector with min and max
vector<vector<int> > nWindow(const int nMolMin, const int nMolMax, const double nExp, const int nWindow, const int nOverlap);
vector<vector<double> > nWindowGrowth(const double mMin, const double mMax, const double grow, const int nWindow, const double dm, const int overlap);
vector<vector<double> > nWindowGrowth(const double mMin, const double mMax, const double grow, const int nWindow, const double dm);

/// given a file name in char*, return basename in string
std::string myTrim(const char* specialchr, const char* fileName);
std::string myTrim(const char* specialchr, string fileName);

/// given a vector of data, return vector of index values for the local maxima
//   tolerance is integer number of neighboring points that must be lower than local max (not counting boundary)
template <class T>
vector<int> findLocalMaxima(const vector<T> data,
  const int itol    //!< tolerance
  ) {
  vector<int> max;

  // process through data
  for (int i = 0; i < int(data.size()); ++i) {

    // determine window range from tolerance, and fixing first and last elements
    int lower = i - itol, upper = i + itol;
    if (lower < 0) lower = 0;
    if (upper >= int(data.size())) upper = int(data.size()) - 1;

    const T windowMax = *std::max_element(data.begin() + lower, data.begin() + upper + 1);
    if (windowMax == data[i]) max.push_back(i);
    //cout << "i " << i << " d " << data[i] << " lower " << lower << " upper " << upper << " wm " << windowMax << " max.size " << max.size() << endl;
  }
  return max;
};

/// given a deque of data, return vector of index values for the local maxima
//   tolerance is integer number of neighboring points that must be lower than local max (not counting boundary)
template <class T>
vector<int> findLocalMaxima(const std::deque<T> data,
  const int itol    //!< tolerance
  ) {
  vector<int> max;

  // process through data
  for (int i = 0; i < int(data.size()); ++i) {

    // determine window range from tolerance, and fixing first and last elements
    int lower = i - itol, upper = i + itol;
    if (lower < 0) lower = 0;
    if (upper >= int(data.size())) upper = int(data.size()) - 1;

    const T windowMax = *std::max_element(data.begin() + lower, data.begin() + upper + 1);
    if (windowMax == data[i]) max.push_back(i);
    //cout << "i " << i << " d " << data[i] << " lower " << lower << " upper " << upper << " wm " << windowMax << " max.size " << max.size() << endl;
  }
  return max;
};

///// given a container (e.g., vector, deque) of data, return vector of index values for the local maxima
////   tolerance is integer number of neighboring points that must be lower than local max (not counting boundary)
//template <template <typename, typename> class Container,
//          typename Value,
//          typename Allocator=std::allocator<Value> >
//vector<int> findLocalMaxima(const Container<Value, Allocator> data,    //!< function
//  const int itol    //!< tolerance
//  ) {
//  vector<int> max;
//
//  // process through data
//  for (int i = 0; i < int(data.size()); ++i) {
//
//    // determine window range from tolerance, and fixing first and last elements
//    int lower = i - itol, upper = i + itol;
//    if (lower < 0) lower = 0;
//    if (upper >= int(data.size())) upper = int(data.size()) - 1;
//
//    const Value windowMax = *std::max_element(data.begin() + lower, data.begin() + upper + 1);
//    if (windowMax == data[i]) max.push_back(i);
//    //cout << "i " << i << " d " << data[i] << " lower " << lower << " upper " << upper << " wm " << windowMax << " max.size " << max.size() << endl;
//  }
//  return max;
//};

// inverse of find local maxima
template <class T>
vector<int> findLocalMinima(const vector<T> data,
  const int itol    //!< tolerance
  ) {
  vector<T> negData = data;
  for (typename vector<T>::iterator it = negData.begin(); it != negData.end(); ++it) *it *= -1;
  return findLocalMaxima(negData, itol);
};

//// inverse of find local maxima
//template <template <typename, typename> class Container,
//          typename Value,
//          typename Allocator=std::allocator<Value> >
//vector<int> findLocalMinima(const Container<Value, Allocator> data,    //!< function
//  const int itol    //!< tolerance
//  ) {
//  Container<Value, Allocator> negData = data;
//  for (typename Container<Value, Allocator>::iterator it = negData.begin(); it != negData.end(); ++it) *it *= -1;
//  return findLocalMaxima(negData, itol);
//};

template<class T>
T mySq(const T x) { return x*x; };

// function to output error messages
void myOut(const char* messageType, std::ostream& message, const std::string className, const int verbose);
void myOutF(const char* messageType, std::ostream& message);

// function to truncate to the nearest factor
template<class T>
T myTrunc(const T x, const T fac) { return x - (x%fac); };

// convert vector of shared pointers to raw pointers
template<class T>
vector<T*> shrPtr2Raw(vector<shared_ptr<T> > shrPtr) {
  vector<T*> raw;
  for (int i = 0; i < int(shrPtr.size()); ++i) {
    raw.push_back(shrPtr[i].get());
  }
  return raw;
};

// return string reported after first appearance of string in file
string fstos(const char* searchString, const char* fileName);
double fstod(const char* searchString, const char* fileName);
int    fstoi(const char* searchString, const char* fileName);
long long fstoll(const char* searchString, const char* fileName);
unsigned long long fstoull(const char* searchString, const char* fileName);

// compute eigenvalues and eigenvectors of a 3x3 real symmetric matrix
int jacobi(vector<vector<double> > matrix, vector<double> &evalues, vector<vector<double> > &evectors);

// preform a single Jacobi rotation
void rotateJacobi(vector<vector<double> > &matrix, const int i, const int j, const int k, const int l, const double s, const double tau);

// compute the squared difference between two vectors
template<class T>
T sqDiff(const vector<T> &x, const vector<T> &y) {
  T sum = 0;
  if (x.size() != y.size()) myOutF("error", std::ostringstream().flush() << "x and y should be same size");
  for (int i = 0; i < int(x.size()); ++i) {
    sum += (x[i] - y[i])*(x[i] - y[i]);
  }
  return sum;
};

/// periodic boundary conditions for orientation angle in 2d
template<class T>
T pbc2d(const T theta) { return -int(theta/2/PI)*2*PI;}; //!< pbc for 2d theta

template<class T>
T pbcangle(const T theta, const T p) { return theta-int(theta/p)*p;}; //!< pbc for 2d theta

/// arc cosine safe from rounding error
template<class T>
T myAcos(T x) { if (x < -1.) x = -1.; if (x > 1.) x = 1; return acos(x); };

/// returns a unit vector orthogonal to given vector
vector<double> orthogonalVec(const vector<double> x);

/// returns rotation matrix for arbitrary rotation axis, u, by an angle theta
vector<vector<double> > rotMatAxisAngle(vector<double> u, const double theta);

/// returns vector rotated by angle theta about axis u
vector<double> rotateVecByAxisAngle(vector<double> x, vector<double> axis, const double theta);

/// draw a superquadric
//  |x/A|^r + |y/B|^s + |z/C|^t <= 1
vector<vector<double> > superQuadric3D(const double A, const double r, const int n, const double umid, const double vmid, const double du, const double dv);
vector<vector<double> > superQuadric3D(const double A, const double r, const int n);
double cosfSuperQuad(const double w, const double m);
double sinfSuperQuad(const double w, const double m);

/// minimum distance between two objects described by grid points
void minDistGrid(const vector<vector<double> > &x1, const vector<vector<double> > &x2, double &xmin, int &x1min, int &x2min);

/// solves quadratic equation, given ax^2+bx+c=0, return x1 and x2. Error if no real solution
template<class T>
T quadraticEqReal(const T a, const T b, const T c, T &x1, T &x2) {
  const T discriminant = b*b-4*a*c;
  if (a*b*c == 0) {
    myOutF("error", std::ostringstream().flush() << "zero coefficient for quadratic equation: a(" << a << ")x^2 + b(" << b << ")x + c(" << c << " ) = 0. Discriminant = " << discriminant << ". Or a coefficient is zero");
  }
  if (discriminant < 0) {
    //myOutF("error", std::ostringstream().flush() << "imaginary roots for quadratic equation: a(" << a << ")x^2 + b(" << b << ")x + c(" << c << " ) = 0. Discriminant = " << discriminant << ". Or a coefficient is zero");
    //    cout << "imaginary roots for quadratic equation: a(" << a << ")x^2 + b(" << b << ")x + c(" << c << " ) = 0. Discriminant = " << discriminant << ". Or a coefficient is zero" << endl;
  } else {
    x1 = (-b+sqrt(discriminant))/2/a;
    x2 = (-b-sqrt(discriminant))/2/a;
  }
  return discriminant;
};

// find value in list, and return true if found, with index of list
template<class T>
bool findInList(const T value, const vector<T> &list, int &index) {
  bool in = false;
  index = -1;
  for (int i = 0; i < int(list.size()); ++i) {
    if (list[i] == value) {
      in = true;
      index = i;
    }
  }
  return in;
};

// find value in list, and return true if found
template<class T>
bool findInList(const T value, const vector<T> &list) { int index; return findInList(value, list, index); };

// convert cartesian coordinates vector to spherical coordinates
vector<double> cartesian2spherical(vector<double> rCart);

// returns the spherical harmonics (l=6) of a cartesian vector r
vector<std::complex<double> > cart2sphereHarm6(vector<double> rCart);

// returns sum of vector of complex numbers multiplied by its conjugate
double complexVec2norm(vector<std::complex<double> > compVec);

// round double to nearest integer
int myRound(double x);

// obtain a rotation matrix from euler angles
vector<vector<double> > Euler2RotMat(const vector<double> euler);

// obtain euler angles from a rotation matrix
vector<vector<double> > RotMat2Euler(const vector<vector<double> > rotMat);

// determinant of a 3x3 matrix
double det3by3(const vector<vector<double> > mat);
double det2by2(const vector<vector<double> > mat);

// cofactor of a 3x3 matrix
vector<vector<double> > cofactor3by3(const vector<vector<double> > mat);

// transpose
vector<vector<double> > transpose(const vector<vector<double> > mat);

// inverse
vector<vector<double> > inv3by3(const vector<vector<double> > mat);

// minimum element in a multidimensional vector
template<class T>
T myMinElement(vector<vector<vector<vector<vector<vector<T> > > > > > vec) {
  vector<T> mins5;
  for (int t = 0; t < int(vec.size()); ++t) {
    vector<T> mins4;
    for (int l = 0; l < int(vec[0].size()); ++l) {
      vector<T> mins3;
      for (int k = 0; k < int(vec[0][0].size()); ++k) {
        vector<T> mins2;
        for (int j = 0; j < int(vec[0][0][0].size()); ++j) {
          vector<T> mins;
          for (int r = 0; r < int(vec[0][0][0][0].size()); ++r) {
            mins.push_back(*std::min_element(vec[t][l][k][j][r].begin(), vec[t][l][k][j][r].begin()+int(vec[t][l][k][j][r].size())));
          }
          mins2.push_back(*std::min_element(mins.begin(), mins.begin()+int(mins.size())));
        }
        mins3.push_back(*std::min_element(mins2.begin(), mins2.begin()+int(mins2.size())));
      }
      mins4.push_back(*std::min_element(mins3.begin(), mins3.begin()+int(mins3.size())));
    }
    mins5.push_back(*std::min_element(mins4.begin(), mins4.begin()+int(mins4.size())));
  }
  return *std::min_element(mins5.begin(), mins5.begin()+int(mins5.size()));
};
template<class T>
T myMaxElement(vector<vector<vector<vector<vector<vector<T> > > > > > vec) {
  vector<T> maxs5;
  for (int t = 0; t < int(vec.size()); ++t) {
    vector<T> maxs4;
    for (int l = 0; l < int(vec[0].size()); ++l) {
      vector<T> maxs3;
      for (int k = 0; k < int(vec[0][0].size()); ++k) {
        vector<T> maxs2;
        for (int j = 0; j < int(vec[0][0][0].size()); ++j) {
          vector<T> maxs;
          for (int r = 0; r < int(vec[0][0][0][0].size()); ++r) {
            maxs.push_back(*std::max_element(vec[t][l][k][j][r].begin(), vec[t][l][k][j][r].begin()+int(vec[t][l][k][j][r].size())));
          }
          maxs2.push_back(*std::max_element(maxs.begin(), maxs.begin()+int(maxs.size())));
        }
        maxs3.push_back(*std::max_element(maxs2.begin(), maxs2.begin()+int(maxs2.size())));
      }
      maxs4.push_back(*std::max_element(maxs3.begin(), maxs3.begin()+int(maxs3.size())));
    }
    maxs5.push_back(*std::max_element(maxs4.begin(), maxs4.begin()+int(maxs4.size())));
  }
  return *std::max_element(maxs5.begin(), maxs5.begin()+int(maxs5.size()));
};

// minimum element in a multidimensional vector
template<class T>
T myMinElement(vector<vector<vector<vector<vector<T> > > > > vec) {
  vector<T> mins4;
  for (int t = 0; t < int(vec.size()); ++t) {
    vector<T> mins3;
    for (int l = 0; l < int(vec[0].size()); ++l) {
      vector<T> mins2;
      for (int k = 0; k < int(vec[0][0].size()); ++k) {
        vector<T> mins;
        for (int j = 0; j < int(vec[0][0][0].size()); ++j) {
          mins.push_back(*std::min_element(vec[t][l][k][j].begin(), vec[t][l][k][j].begin()+int(vec[t][l][k][j].size())));
        }
        mins2.push_back(*std::min_element(mins.begin(), mins.begin()+int(mins.size())));
      }
      mins3.push_back(*std::min_element(mins2.begin(), mins2.begin()+int(mins2.size())));
    }
    mins4.push_back(*std::min_element(mins3.begin(), mins3.begin()+int(mins3.size())));
  }
  return *std::min_element(mins4.begin(), mins4.begin()+int(mins4.size()));
};
template<class T>
T myMaxElement(vector<vector<vector<vector<vector<T> > > > > vec) {
  vector<T> maxs4;
  for (int t = 0; t < int(vec.size()); ++t) {
    vector<T> maxs3;
    for (int l = 0; l < int(vec[0].size()); ++l) {
      vector<T> maxs2;
      for (int k = 0; k < int(vec[0][0].size()); ++k) {
        vector<T> maxs;
        for (int j = 0; j < int(vec[0][0][0].size()); ++j) {
          maxs.push_back(*std::max_element(vec[t][l][k][j].begin(), vec[t][l][k][j].begin()+int(vec[t][l][k][j].size())));
        }
        maxs2.push_back(*std::max_element(maxs.begin(), maxs.begin()+int(maxs.size())));
      }
      maxs3.push_back(*std::max_element(maxs2.begin(), maxs2.begin()+int(maxs2.size())));
    }
    maxs4.push_back(*std::max_element(maxs3.begin(), maxs3.begin()+int(maxs3.size())));
  }
  return *std::max_element(maxs4.begin(), maxs4.begin()+int(maxs4.size()));
};

// minimum element in a multidimensional vector
template<class T>
T myMinElement(vector<vector<vector<vector<T> > > > vec) {
  vector<T> mins3;
  for (int l = 0; l < int(vec.size()); ++l) {
    vector<T> mins2;
    for (int k = 0; k < int(vec[0].size()); ++k) {
      vector<T> mins;
      for (int j = 0; j < int(vec[0][0].size()); ++j) {
        mins.push_back(*std::min_element(vec[l][k][j].begin(), vec[l][k][j].begin()+int(vec[l][k][j].size())));
      }
      mins2.push_back(*std::min_element(mins.begin(), mins.begin()+int(mins.size())));
    }
    mins3.push_back(*std::min_element(mins2.begin(), mins2.begin()+int(mins2.size())));
  }
  return *std::min_element(mins3.begin(), mins3.begin()+int(mins3.size()));
};
template<class T>
T myMaxElement(vector<vector<vector<vector<T> > > > vec) {
  vector<T> maxs3;
  for (int l = 0; l < int(vec.size()); ++l) {
    vector<T> maxs2;
    for (int k = 0; k < int(vec[0].size()); ++k) {
      vector<T> maxs;
      for (int j = 0; j < int(vec[0][0].size()); ++j) {
        maxs.push_back(*std::max_element(vec[l][k][j].begin(), vec[l][k][j].begin()+int(vec[l][k][j].size())));
      }
      maxs2.push_back(*std::max_element(maxs.begin(), maxs.begin()+int(maxs.size())));
    }
    maxs3.push_back(*std::max_element(maxs2.begin(), maxs2.begin()+int(maxs2.size())));
  }
  return *std::max_element(maxs3.begin(), maxs3.begin()+int(maxs3.size()));
};
// minimum element in a multidimensional vector
template<class T>
T myMinElement(vector<vector<vector<T> > > vec) {
  vector<T> mins2;
  for (int k = 0; k < int(vec.size()); ++k) {
    vector<T> mins;
    for (int j = 0; j < int(vec[0].size()); ++j) {
      mins.push_back(*std::min_element(vec[k][j].begin(), vec[k][j].begin()+int(vec[k][j].size())));
    }
    mins2.push_back(*std::min_element(mins.begin(), mins.begin()+int(mins.size())));
  }
  return *std::min_element(mins2.begin(), mins2.begin()+int(mins2.size()));
};
template<class T>
T myMaxElement(vector<vector<vector<T> > > vec) {
  vector<T> maxs2;
  for (int k = 0; k < int(vec.size()); ++k) {
    vector<T> maxs;
    for (int j = 0; j < int(vec[0].size()); ++j) {
      maxs.push_back(*std::max_element(vec[k][j].begin(), vec[k][j].begin()+int(vec[k][j].size())));
    }
    maxs2.push_back(*std::max_element(maxs.begin(), maxs.begin()+int(maxs.size())));
  }
  return *std::max_element(maxs2.begin(), maxs2.begin()+int(maxs2.size()));
};

/// convert quaternion vector to euler angles (uses quat2rot -> rot2euler)
vector<double> quat2euler(vector<double> quat);

/// convert euler angles to quaternions
// NOTE: This doesn't work!!
vector<double> rot2quat(const vector<vector<double> > &rot);

/// periodic boundary conditions for euler angles
void eulerPBC(vector<double> &euler);

/// output matrix
template<class T>
void myCout(vector<T> data) {
  for (unsigned int i = 0; i < data.size(); ++i) {
    cout << data[i] << " ";
  }
  cout << endl;
};
template<class T>
void myCout(vector<vector<T> > data) {
  for (unsigned int i = 0; i < data.size(); ++i) {
    for (unsigned int j = 0; j < data.size(); ++j) {
      cout << data[i][j] << " ";
    }
    cout << endl;
  }
};

/// my acos includes some double precision checks
double myACOS(double x);

/// trapezoid rule integration
template<class T>
double myIntegratorTrapezoid(vector<T> data) {
  double accumulator = 0;
  for (unsigned int i = 0; i < data.size()-1; ++i) {
    accumulator += 0.5*(data[i]+data[i+1]);
  }
  return accumulator;
}

/// resize vectors
template<class T>
void myResize(const int xs, const int ys, vector<vector<T> > *vec) {
  vec->resize(xs);
  for (unsigned int i = 0; i < vec->size(); ++i) {
    (*vec)[i].resize(ys);
  }
}

/// outer product of two vectors returns matrix
template<class T>
vector<vector<T> > outerProd(const vector<T> &u, const vector<T> &v) {
  vector<vector<T> > res(u.size(), vector<T>(v.size()));
  for (unsigned int i = 0; i < u.size(); ++i) {
    for (unsigned int j = 0; j < v.size(); ++j) {
      res[i][j] = u[i]*v[j];
    }
  }
  return res;
}

#endif  //FUNCTIONS_H_
