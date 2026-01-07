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

  // Insert a nested table.
  void insert(const int bin, const RecursiveTable1D& nested);

  double linear_interpolation(const double value0) const override;
  double forward_difference_interpolation(const double value0) const override;

  void serialize(std::ostream& ostr) const;
  explicit RecursiveTable1D(std::istream& istr);
  virtual ~RecursiveTable1D();

 private:
  std::vector<std::unique_ptr<RecursiveTable1D> > nested_;
};

}  // namespace feasst

#endif  // FEASST_MATH_RECURSIVE_TABLE_H_
