/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of this agreement (see LICENSE.txt) and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

/**
 * \file
 *
 * \brief data table with interpolation
 */
#ifndef TABLE_H_
#define TABLE_H_

#include <string>
#include <vector>
#ifdef GSL_
  #include <stdlib.h>
  #include <stdio.h>
  #include <math.h>
  #include <gsl/gsl_errno.h>
  #include <gsl/gsl_spline.h>
#endif  // GSL_
#include "./base.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class Table {
 public:
  Table();
  explicit Table(const char* fileName);
  virtual ~Table();
  virtual Table* clone() const { Table* p = new Table(*this); return p; }

  /// linear interpolation
  double interpolate(const double val0);
  double interpolate(const double val0, double *deriv);
  double interpolate(const double val0, const double val1);
  double interpolate(const double val0, const double val1, const double val2);
  double interpolate(const double val0, const double val1, const double val2,
                     const double val3);
  double interpolate(const double val0, const double val1, const double val2,
                     const double val3, const double val4);
  double interpolate(const double val0, const double val1, const double val2,
                     const double val3, const double val4, const double val5);

  /// compute minimum value in table
  double compute_min() const;
  double compute_max() const;

  /// compute minimium values as a function of one table dimension, dim
  void compute_min_compress1d(const int dim);

  // return minimum value by interpolation of the first dimension
  double min(const double val0) { return interpolate(val0); }

  /// print table in hdf5 format
  void printHDF5(const char* fileName);

  /// set the interpolator
  void setInterpolator(const char* name);

  /// solve for the spline derivatives, assuming condition at end points
  void solveSpline(const char* endCondition);

  /// for a given bin, return the abscissae ("x")
  double bin2abs(const int bin);

  // functions for read-only access of private data-members
  double min() const { return min_; }
  vector<vector<double> > tablim() const { return tablim_; }
  vector<vector<vector<double> > > tab3() const { return tab3_; }
  vector<vector<vector<vector<double> > > > tab4() const { return tab4_; }
  string interpolator() const { return interpolator_; }

 protected:
  int tabDims_;             //!< number of dimensions in table
  string tabType_;          //!< type of table
  string interpolator_;     //!< type of interpolation
  vector<double> tab1_;     //!< table 1D
  vector<vector<double> > tab2_;   //!< table 2D
  vector<vector<vector<double> > > tab3_;   //!< table 3D
  vector<vector<vector<vector<double> > > > tab4_;   //!< table 4D
  vector<vector<vector<vector<vector<double> > > > > tab5_;   //!< table 5D
  vector<vector<vector<vector<vector<vector<double> > > > > > tab6_;
  vector<vector<double> > tablim_;   //!< table limits
  double min_;              //!< minimum value
  double d0_, d1_, d2_, d3_, d4_, d5_;

  // spline
  vector<double> cspline_;   //!< spline coefficients

  /// accessor function with indices [1,n] to look like FORTRAN
  double c_(const int coeff, const int index) const {
    return cspline_[4*(index-1)+(coeff-1)];
  }
  void cset_(const int coeff, const int index, const double value) {
    cspline_[4*(index-1)+(coeff-1)] = value;
  }

  #ifdef GSL_
    gsl_interp_accel *acc;
    gsl_spline *spline;
  #endif  // GSL_

  // defaults in constructor
  void defaultConstruction_();
};

class erftable {
 public:
  erftable() {}
  ~erftable() {}
  void init(const double alpha, const double rCut);
  double eval(const double x) const;
 private:
  vector<double> vtab_;
  int n_;
  double ds_;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TABLE_H_

