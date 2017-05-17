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
#ifdef _OPENMP
  #include <omp.h>
#endif  // _OPENMP
#include "./custom_exception.h"
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::shared_ptr;
using std::make_shared;

/// If the assertion condition is not true, throw exception with message.
# define ASSERT(condition, message) \
if (! (condition)) { \
  std::stringstream err_msg; \
  err_msg << "Assertion `" #condition "` failed in " << __FILE__ \
            << " line " << __LINE__ << ": " << message; \
  CustomException c(err_msg); \
}

/// If the warning condition is true, send message to standard output.
# define WARN(condition, message) \
if (condition) { \
  std::cout << "Warning `" #condition "` in " << __FILE__ \
            << " line " << __LINE__ << ": " << message << std::endl; \
}

/// Send message to standard output.
# define NOTE(message) \
std::cout << "Note in " << __FILE__ \
          << " line " << __LINE__ << ": " << message << std::endl; \
//throw c;

namespace feasst {

/// Return the magnitude of the first argument, with the sign of the second.
double sign(const double a, const double b);

/// Function to fill 2D vector with input variable.
void fill(const double input, vector<vector<double> > &x);
void fill(const int input, vector<vector<int> > &x);

/// Function to fill 3D vector with input variable.
template<class T>
void fill(const T input, vector<vector<vector<T> > > &x) {
  for (typename vector<vector<vector<T> > >::iterator iter = x.begin();
       iter != x.end(); ++iter) {
    for (typename vector<vector<T> >::iterator iter2 = (*iter).begin();
         iter2 != (*iter).end(); ++iter2) {
      std::fill((*iter2).begin(), (*iter2).end(), input);
    }
  }
};

/// \return the number of lines in a file
int numLines(const string fileName);

/// \return vector which is product of matrix and vector e.g., A[][] x[] = b[]
vector<double> matVecMul(const vector<vector<double> > &a, 
  const vector<double> &x);

/// \return vector (inner) scalar dot product, a[] . b[] = scalar
double vecDotProd(const vector<double> &a, const vector<double> &b);

/// normalize vector to unit size
void normalizeVec(vector<double> *x);

/// \return vector cross product
vector<double> crossProd(const vector<double> &a, const vector<double> &b);

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

/// initialize random number generator based on date
void ranInitByDate();

/// initialize random number generator to seed value for reproducibility
void ranInitForRepro(const int seed = 1346867550);

/**
 * \return rotation matrix from quaternions
 * Franz J. Vesely, J. Comput. Phys., 47, 291-296 (1982)
 */
vector<vector<double> > quat2rot(vector<double> quaternion);

/// \return 2d rotation matrix from angle
vector<vector<double> > theta2rot(double theta);

/// \return product of all elements of a vector
template<class T>
T product(const vector<T> &vec) {
  T prod = 1;
  for (int i = 0; i < int(vec.size()); ++i) {
    prod *= vec[i];
  }
  return prod;
};

/**
 *  \return vector x with atomic positions of a reference SPC/E water molecule
 *  vector of atoms, oxygen first
 */
vector<vector<double> > vecSPCE();

/// \return product of two matrices
vector<vector<double> > matMul(const vector<vector<double> > &a,
  const vector<vector<double> > &b);

/// \return volume of spherical shell
double volShell(const double rabove,     //!< upper limit of shell
                const double rbelow,     //!< lower limit of shell
                const int dim = 3        //!< spatial dimensions
);

/// \return average of vector
template<class T>
double vecAverage(const vector<T> &x) {
  return std::accumulate(x.begin(), x.end(), 0.) / int(x.size());
};

/// \return true if file exists
bool fileExists(const char* fileName);

/// \return true if file exists
bool fileExists(std::ifstream& file);

/// renames file by appending with ".bak"
void fileBackUp(const char* fileName);

/// renames file by appending with ".bak"
void fileBackUp(const std::string fileName);

/// skip all lines beginning with character in file
void skipCharsInFile(const char comment, std::ifstream &file);

/// function to skip all lines until reaching a certain line
void readUntil(const char* searchString, std::ifstream &file);

/**
 * given range of n, exponent, the number of windows, and overlap of windows, return vector with min and max
 */
vector<vector<int> > nWindow(const int nMolMin,   //!< minimum n
  const int nMolMax,   //!< maximum n
  const double nExp,   //!< exponent of distribution
  const int nWindow,   //!< number of windows
  const int nOverlap   //!< overlap of n between windows
);

/**
 * \return the range of windows in parallelization of an order parameter
 *  the order parameter of the windows is in the range [mMin, mMax]
 *  the number of windows is nWindow
 *  the bin width of the order parameter is dm
 *
 *  the windows vary in sized based on the grow parameter.
 *  if grow == 0, the windows are roughly the same size
 *  if grow <  0, the windows decrease in size
 *  if grow >  0, the windows increase in size
 *  Note: grow != 0 is not working correctly
 *
 *  For example, if grow > 0, n=3, the windows may look like this
 *
 *  |---------------------| total window
 *  |--|------|-----------| three windows
 *
 */
vector<vector<double> > nWindowGrowth(const double mMin,   //!< minimum m
  const double mMax,     //!< maximum m
  const double grow,     //!< percentage to grow each window
  const int nWindow,     //!< number of windows
  const double dm,       //!< width of m
  const int overlap = 0  //!< overlap+1 overlapping windows
);

/// \return fileName with all characters up to the last specialchr removed
std::string trim(const char* specialchr, const char* fileName);

/// \return fileName with all characters up to the last specialchr removed
std::string trim(const char* specialchr, string fileName);

/// \return vector of index values for the local maxima.
template <class T>
vector<int> findLocalMaxima(const vector<T> data,
  const int itol  /** integer number of neighboring points that must be lower
    than local max (not counting boundary) */
  ) {
  vector<int> max;

  // process through data
  for (int i = 0; i < int(data.size()); ++i) {

    // determine window range from tolerance, and fixing first and last elements
    int lower = i - itol, upper = i + itol;
    if (lower < 0) lower = 0;
    if (upper >= int(data.size())) upper = int(data.size()) - 1;

    const T windowMax = *std::max_element(data.begin() + lower,
      data.begin() + upper + 1);
    if (windowMax == data[i]) max.push_back(i);
  }
  return max;
};

/**
 *  \return vector of index values for the local maxima.
 * HWH NOTE: This is copy and pasted from above with vector, but 
 * implementation with multiple templates leads to errors in swig. 
 * See commented implementation of findLocalMinimum below.
 */
template <class T>
vector<int> findLocalMaxima(const std::deque<T> data,
  const int itol  /** integer number of neighboring points that must be lower
    than local max (not counting boundary) */
  ) {
  vector<int> max;

  // process through data
  for (int i = 0; i < int(data.size()); ++i) {

    // determine window range from tolerance, and fixing first and last elements
    int lower = i - itol, upper = i + itol;
    if (lower < 0) lower = 0;
    if (upper >= int(data.size())) upper = int(data.size()) - 1;

    const T windowMax = *std::max_element(data.begin() + lower,
      data.begin() + upper + 1);
    if (windowMax == data[i]) max.push_back(i);
  }
  return max;
};

/// \return inverse of find local maxima. Simply multiply the data by -1
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

/// \return vector of raw pointers from vector of shared pointers
template<class T>
vector<T*> shrPtr2Raw(vector<shared_ptr<T> > shrPtr) {
  vector<T*> raw;
  for (int i = 0; i < int(shrPtr.size()); ++i) {
    raw.push_back(shrPtr[i].get());
  }
  return raw;
};

/// \return string reported after first appearance of searchString in file
string fstos(const char* searchString, const char* fileName);

/// \return double reported after first appearance of searchString in file
double fstod(const char* searchString, const char* fileName);

/// \return integer reported after first appearance of searchString in file
int fstoi(const char* searchString, const char* fileName);

/// \return long long reported after first appearance of searchString in file
long long fstoll(const char* searchString, const char* fileName);

/**
 * \return unsigned long long reported after first appearance of searchString
 *  in file
 */
unsigned long long fstoull(const char* searchString, const char* fileName);

/**
 * compute eigenvalues and eigenvectors based on Jacobi rotations. Adapted 
 * from LAMMPS, which was adapted from Numerical Recipes jacobi() function
 *  \return error code (1 if failed) 
 */
int jacobi(vector<vector<double> > matrix,  //!< 3x3 real symmetric matrix
  vector<double> &evalues,            //!< computed eigen values
  vector<vector<double> > &evectors   //!< computed eigen vectors
);

/**
 * preform a single Jacobi rotation. Adapted from LAMMPS, via NR.
 */
void rotateJacobi(vector<vector<double> > &matrix, const int i, const int j,
  const int k, const int l, const double s, const double tau);

// \return the squared difference between two vectors
template<class T>
T sqDiff(const vector<T> &x, const vector<T> &y) {
  ASSERT(x.size() == y.size(), "x and y should be same size");
  T sum = 0;
  for (int i = 0; i < int(x.size()); ++i) {
    sum += (x[i] - y[i])*(x[i] - y[i]);
  }
  return sum;
};

/// \return periodic boundary conditions for orientation angle in 2d
template<class T>
T pbc2d(const T theta) { return -int(theta/2./PI)*2.*PI;};

/// \return arc cosine which is safe from rounding error
template<class T>
T arccos(T x) { if (x < -1.) x = -1.; if (x > 1.) x = 1; return acos(x); };

/// \return unit vector orthogonal to given vector
vector<double> orthogonalVec(const vector<double> &x);

/// \return rotation matrix for axis angle rotation
vector<vector<double> > rotMatAxisAngle(
  vector<double> axis, const double theta);

/// \return vector rotated by angle theta about axis
vector<double> rotateVecByAxisAngle(vector<double> x, vector<double> axis,
  const double theta);

/**
 * Solves quadratic equation, given ax^2+bx+c=0.
 * \return discriminant b^2-4ac
 */
template<class T>
T quadraticEqReal(const T a, const T b, const T c,
  T &x1,  //!< first root of the quadratic equation
  T &x2   //!< second root of the quadratic equation
  ) {
  const T discriminant = b*b-4*a*c;
  ASSERT(a*b*c != 0, "zero coefficient for quadratic equation: a(" << a 
    << ")x^2 + b(" << b << ")x + c(" << c << " ) = 0. Discriminant = " 
    << discriminant << ". Or a coefficient is zero");
  if (discriminant >= 0) {
    x1 = (-b+sqrt(discriminant))/2/a;
    x2 = (-b-sqrt(discriminant))/2/a;
  }
  return discriminant;
};

/// \return if value is found in list
template<class T>
bool findInList(const T value, const vector<T> &list,
  int &index  //!< last index in list where value was found
  ) {
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

/// \return if value is found in list
template<class T>
bool findInList(const T value, const vector<T> &list) { 
  int index;
  return findInList(value, list, index); 
};

/// \return spherical coordinates given cartesian coordinates
vector<double> cartesian2spherical(vector<double> rCart);

/// \return spherical harmonics (l=6) of a cartesian vector r
vector<std::complex<double> > cart2sphereHarm6(vector<double> rCart);

/// \return sum of vector of complex numbers multiplied by its conjugate
double complexVec2norm(vector<std::complex<double> > compVec);

/// \return rounded double to nearest integer
int round(double x);

/**
 * \return rotation matrix from euler angles
 *
 * http://mathworld.wolfram.com/EulerAngles.html
 *
 * Euler angles (phi, theta, psi) are defined by rotation phi about z-axis,
 *  rotation theta about new x-xis, then rotation psi about new z-axis
 *  phi [-pi,pi], theta [0, pi], psi [-pi, pi]
 */
vector<vector<double> > Euler2RotMat(const vector<double> euler);

/// \return euler angles from rotation matrix (same convention as Euler2RotMat)
vector<vector<double> > RotMat2Euler(const vector<vector<double> > rotMat);

/// \return determinant of a 3x3 matrix
double det3by3(const vector<vector<double> > matrix);

/// \return determinant of a 2x2 matrix
double det2by2(const vector<vector<double> > matrix);

/// \return cofactor of a 3x3 matrix
vector<vector<double> > cofactor3by3(const vector<vector<double> > matrix);

/// \return transpose of a matrix
template<class T>
vector<vector<T> > transpose(const vector<vector<T> > matrix) {
  vector<vector<T> > transpose(int(matrix[0].size()), vector<T>(
    int(matrix.size())));
  for (int idim = 0; idim < int(matrix.size()); ++idim) {
    for (int jdim = 0; jdim < int(matrix[idim].size()); ++jdim) {
      transpose[jdim][idim] = matrix[idim][jdim];
    }
  }
  return transpose;
}


/// \return inverse of a 3x3 matrix
vector<vector<double> > inv3by3(const vector<vector<double> > matrix);

/// \return minimum element in a multidimensional vector
template<class T>
T minElement(vector<vector<vector<vector<vector<vector<T> > > > > > vec) {
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
            mins.push_back(*std::min_element(vec[t][l][k][j][r].begin(),
              vec[t][l][k][j][r].begin()+int(vec[t][l][k][j][r].size())));
          }
          mins2.push_back(*std::min_element(mins.begin(),
            mins.begin()+int(mins.size())));
        }
        mins3.push_back(*std::min_element(mins2.begin(),
          mins2.begin()+int(mins2.size())));
      }
      mins4.push_back(*std::min_element(mins3.begin(),
        mins3.begin()+int(mins3.size())));
    }
    mins5.push_back(*std::min_element(mins4.begin(),
      mins4.begin()+int(mins4.size())));
  }
  return *std::min_element(mins5.begin(), mins5.begin()+int(mins5.size()));
};

/// \return maximum element in a multidimensional vector
template<class T>
T maxElement(vector<vector<vector<vector<vector<vector<T> > > > > > vec) {
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
            maxs.push_back(*std::max_element(vec[t][l][k][j][r].begin(),
              vec[t][l][k][j][r].begin()+int(vec[t][l][k][j][r].size())));
          }
          maxs2.push_back(*std::max_element(maxs.begin(),
            maxs.begin()+int(maxs.size())));
        }
        maxs3.push_back(*std::max_element(maxs2.begin(),
          maxs2.begin()+int(maxs2.size())));
      }
      maxs4.push_back(*std::max_element(maxs3.begin(),
        maxs3.begin()+int(maxs3.size())));
    }
    maxs5.push_back(*std::max_element(maxs4.begin(),
      maxs4.begin()+int(maxs4.size())));
  }
  return *std::max_element(maxs5.begin(), maxs5.begin()+int(maxs5.size()));
};

template<class T>
T minElement(vector<vector<vector<vector<vector<T> > > > > vec) {
  vector<T> mins4;
  for (int t = 0; t < int(vec.size()); ++t) {
    vector<T> mins3;
    for (int l = 0; l < int(vec[0].size()); ++l) {
      vector<T> mins2;
      for (int k = 0; k < int(vec[0][0].size()); ++k) {
        vector<T> mins;
        for (int j = 0; j < int(vec[0][0][0].size()); ++j) {
          mins.push_back(*std::min_element(vec[t][l][k][j].begin(),
            vec[t][l][k][j].begin()+int(vec[t][l][k][j].size())));
        }
        mins2.push_back(*std::min_element(mins.begin(),
          mins.begin()+int(mins.size())));
      }
      mins3.push_back(*std::min_element(mins2.begin(),
        mins2.begin()+int(mins2.size())));
    }
    mins4.push_back(*std::min_element(mins3.begin(),
      mins3.begin()+int(mins3.size())));
  }
  return *std::min_element(mins4.begin(), mins4.begin()+int(mins4.size()));
};
template<class T>
T maxElement(vector<vector<vector<vector<vector<T> > > > > vec) {
  vector<T> maxs4;
  for (int t = 0; t < int(vec.size()); ++t) {
    vector<T> maxs3;
    for (int l = 0; l < int(vec[0].size()); ++l) {
      vector<T> maxs2;
      for (int k = 0; k < int(vec[0][0].size()); ++k) {
        vector<T> maxs;
        for (int j = 0; j < int(vec[0][0][0].size()); ++j) {
          maxs.push_back(*std::max_element(vec[t][l][k][j].begin(),
            vec[t][l][k][j].begin()+int(vec[t][l][k][j].size())));
        }
        maxs2.push_back(*std::max_element(maxs.begin(),
          maxs.begin()+int(maxs.size())));
      }
      maxs3.push_back(*std::max_element(maxs2.begin(),
        maxs2.begin()+int(maxs2.size())));
    }
    maxs4.push_back(*std::max_element(maxs3.begin(),
      maxs3.begin()+int(maxs3.size())));
  }
  return *std::max_element(maxs4.begin(), maxs4.begin()+int(maxs4.size()));
};

template<class T>
T minElement(vector<vector<vector<vector<T> > > > vec) {
  vector<T> mins3;
  for (int l = 0; l < int(vec.size()); ++l) {
    vector<T> mins2;
    for (int k = 0; k < int(vec[0].size()); ++k) {
      vector<T> mins;
      for (int j = 0; j < int(vec[0][0].size()); ++j) {
        mins.push_back(*std::min_element(vec[l][k][j].begin(),
          vec[l][k][j].begin()+int(vec[l][k][j].size())));
      }
      mins2.push_back(*std::min_element(mins.begin(),
        mins.begin()+int(mins.size())));
    }
    mins3.push_back(*std::min_element(mins2.begin(),
      mins2.begin()+int(mins2.size())));
  }
  return *std::min_element(mins3.begin(), mins3.begin()+int(mins3.size()));
};
template<class T>
T maxElement(vector<vector<vector<vector<T> > > > vec) {
  vector<T> maxs3;
  for (int l = 0; l < int(vec.size()); ++l) {
    vector<T> maxs2;
    for (int k = 0; k < int(vec[0].size()); ++k) {
      vector<T> maxs;
      for (int j = 0; j < int(vec[0][0].size()); ++j) {
        maxs.push_back(*std::max_element(vec[l][k][j].begin(),
          vec[l][k][j].begin()+int(vec[l][k][j].size())));
      }
      maxs2.push_back(*std::max_element(maxs.begin(),
        maxs.begin()+int(maxs.size())));
    }
    maxs3.push_back(*std::max_element(maxs2.begin(),
      maxs2.begin()+int(maxs2.size())));
  }
  return *std::max_element(maxs3.begin(), maxs3.begin()+int(maxs3.size()));
};
template<class T>
T minElement(vector<vector<vector<T> > > vec) {
  vector<T> mins2;
  for (int k = 0; k < int(vec.size()); ++k) {
    vector<T> mins;
    for (int j = 0; j < int(vec[0].size()); ++j) {
      mins.push_back(*std::min_element(vec[k][j].begin(),
        vec[k][j].begin()+int(vec[k][j].size())));
    }
    mins2.push_back(*std::min_element(mins.begin(),
      mins.begin()+int(mins.size())));
  }
  return *std::min_element(mins2.begin(), mins2.begin()+int(mins2.size()));
};
template<class T>
T maxElement(vector<vector<vector<T> > > vec) {
  vector<T> maxs2;
  for (int k = 0; k < int(vec.size()); ++k) {
    vector<T> maxs;
    for (int j = 0; j < int(vec[0].size()); ++j) {
      maxs.push_back(*std::max_element(vec[k][j].begin(),
        vec[k][j].begin()+int(vec[k][j].size())));
    }
    maxs2.push_back(*std::max_element(maxs.begin(),
      maxs.begin()+int(maxs.size())));
  }
  return *std::max_element(maxs2.begin(), maxs2.begin()+int(maxs2.size()));
};

/// \return euler angles given quaternions
vector<double> quat2euler(vector<double> quat);

/// resize vectors
template<class T>
void resize(const int xs, const int ys, vector<vector<T> > *vec) {
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

/// \return string representation of vector
template<class T>
std::string vec2str(const vector<T> &vec) {
  std::stringstream ss;
  for (unsigned int i = 0; i < vec.size(); ++i) {
    ss << vec[i] << " ";
  }
  return ss.str();
}
template<class T>
std::string vec2str(const vector<vector<T> > &vec) {
  std::stringstream ss;
  for (unsigned int i = 0; i < vec.size(); ++i) {
    ss << vec2str(vec[i]) << endl;
  }
  return ss.str();
}

} // namespace feasst

#endif  //FUNCTIONS_H_
