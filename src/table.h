/**
 * \file
 *
 * \brief data table with interpolation
 */
#ifndef TABLE_H_
#define TABLE_H_

#include <string>
#include <vector>
#include "./base_all.h"

class Table {
 public:
  Table();
  explicit Table(const char* fileName);
  virtual ~Table() {}
  virtual Table* clone() const { Table* p = new Table(*this); return p; }

  // defaults in constructor
  void defaultConstruction();

  /// linear interpolation
  double interpolate(const double val0);
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

  // functions for read-only access of private data-members
  double min() const { return min_; }
  vector<vector<double> > tablim() const { return tablim_; }
  vector<vector<vector<double> > > tab3() const { return tab3_; }
  vector<vector<vector<vector<double> > > > tab4() const { return tab4_; }

 protected:
  int tabDims_;             //!< number of dimensions in table
  string tabType_;          //!< type of table
  vector<double> tab1_;     //!< table 1D
  vector<vector<double> > tab2_;   //!< table 2D
  vector<vector<vector<double> > > tab3_;   //!< table 3D
  vector<vector<vector<vector<double> > > > tab4_;   //!< table 4D
  vector<vector<vector<vector<vector<double> > > > > tab5_;   //!< table 5D
  vector<vector<vector<vector<vector<vector<double> > > > > > tab6_;
  vector<vector<double> > tablim_;   //!< table limits
  double min_;              //!< minimum value
  double d0_, d1_, d2_, d3_, d4_, d5_;
};

#endif  // TABLE_H_

