#ifndef FEASST_MATH_RECURSIVE_TABLE_H_
#define FEASST_MATH_RECURSIVE_TABLE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "math/include/table.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  One-dimensional nested table.

  The bounds within a nested table are repeated, so it is inefficient if nested
  tables are small.
 */
class RecursiveTable1D : public Table1D {
 public:
  /**
    args:
    - Table1D arguments.
   */
  explicit RecursiveTable1D(argtype args = argtype());
  explicit RecursiveTable1D(argtype *args);

  /// Insert a nested table.
  void insert(const int bin, const RecursiveTable1D& nested);

  /// Return the percentage of table that is nested.
  double percent_nested() const;

  double linear_interpolation(const double value0) const override;
  double forward_difference_interpolation(const double value0) const override;

  void serialize(std::ostream& ostr) const;
  explicit RecursiveTable1D(std::istream& istr);
  virtual ~RecursiveTable1D();

 private:
  std::vector<std::shared_ptr<RecursiveTable1D> > nested_;
  //std::vector<std::unique_ptr<RecursiveTable1D> > nested_;
};

/**
  Two-dimensional nested table.
 */
class RecursiveTable2D : public Table2D {
 public:
  /**
    args:
    - Table2D arguments.
   */
  explicit RecursiveTable2D(argtype args = argtype());
  explicit RecursiveTable2D(argtype *args);

  // Insert a nested table.
  void insert(const int bin0, const int bin1, const RecursiveTable2D& nested);

  double linear_interpolation(const double value0,
                              const double value1) const override;

  /// Return the total number of data points.
  int num_data() const;

  /// Return the percentage of table that is nested.
  double percent_nested() const;

  void serialize(std::ostream& ostr) const;
  explicit RecursiveTable2D(std::istream& istr);
  virtual ~RecursiveTable2D();

 private:
  std::vector<std::vector<std::shared_ptr<RecursiveTable2D> > > nested_;
};

/**
  Three-dimensional nested table.
 */
class RecursiveTable3D : public Table3D {
 public:
  /**
    args:
    - Table3D arguments.
   */
  explicit RecursiveTable3D(argtype args = argtype());
  explicit RecursiveTable3D(argtype *args);

  // Insert a nested table.
  void insert(const int bin0, const int bin1, const int bin2,
    const RecursiveTable3D& nested);

  /// Return the percentage of table that is nested.
  double percent_nested() const;

  double linear_interpolation(const double value0, const double value1,
    const double value2) const override;

  void serialize(std::ostream& ostr) const;
  explicit RecursiveTable3D(std::istream& istr);
  virtual ~RecursiveTable3D();

 private:
  std::vector<std::vector<std::vector<std::shared_ptr<RecursiveTable3D> > > > nested_;
};

/**
  Five-dimensional nested table.
 */
class RecursiveTable5D : public Table5D {
 public:
  /**
    args:
    - Table5D arguments.
   */
  explicit RecursiveTable5D(argtype args = argtype());
  explicit RecursiveTable5D(argtype *args);

  // Insert a nested table.
  void insert(const int bin0, const int bin1, const int bin2, const int bin3,
    const int bin4, const RecursiveTable5D& nested);

  double linear_interpolation(const double value0, const double value1,
    const double value2, const double value3,
    const double value4) const override;

  void serialize(std::ostream& ostr) const;
  explicit RecursiveTable5D(std::istream& istr);
  virtual ~RecursiveTable5D();

 private:
  std::vector<std::vector<std::vector<std::vector<std::vector<std::shared_ptr<RecursiveTable5D> > > > > > nested_;
};

/**
  Six-dimensional nested table.
 */
class RecursiveTable6D : public Table6D {
 public:
  /**
    args:
    - Table6D arguments.
   */
  explicit RecursiveTable6D(argtype args = argtype());
  explicit RecursiveTable6D(argtype *args);

  // Insert a nested table.
  void insert(const int bin0, const int bin1, const int bin2, const int bin3,
    const int bin4, const int bin5, const RecursiveTable6D& nested);

  double linear_interpolation(const double value0, const double value1,
    const double value2, const double value3,
    const double value4, const double value5) const override;

  void serialize(std::ostream& ostr) const;
  explicit RecursiveTable6D(std::istream& istr);
  virtual ~RecursiveTable6D();

 private:
  std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::shared_ptr<RecursiveTable6D> > > > > > > nested_;
};

}  // namespace feasst

#endif  // FEASST_MATH_RECURSIVE_TABLE_H_
