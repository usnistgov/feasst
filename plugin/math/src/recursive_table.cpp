
#include "utils/include/utils.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/serialize_extra.h"
#include "math/include/recursive_table.h"

namespace feasst {

RecursiveTable1D::RecursiveTable1D(argtype * args) : Table1D(args) {
  nested_.resize(num());
}
RecursiveTable1D::RecursiveTable1D(argtype args) : RecursiveTable1D(&args) { feasst_check_all_used(args); }
RecursiveTable1D::~RecursiveTable1D() {}

void RecursiveTable1D::serialize(std::ostream& ostr) const {
  Table1D::serialize(ostr);
  feasst_serialize_version(1864, ostr);
  feasst_serialize(nested_, ostr);
}

RecursiveTable1D::RecursiveTable1D(std::istream& istr) : Table1D(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1864, "version: " << version);
//  feasst_deserialize(&nested_, istr);
   int dim1;
   istr >> dim1;
   nested_.resize(dim1);
   for (int index = 0; index < dim1; ++index) {
     //feasst_deserialize((*vector)[index], istr);
     int existing;
     istr >> existing;
     if (existing != 0) {
       nested_[index] = std::make_shared<RecursiveTable1D>(istr);
     }
   }
}

void RecursiveTable1D::insert(const int bin, const RTable& nested) {
  ASSERT(bin >= 0 && bin < num(), "bin:" << bin << " must be >0 and < num:"
    << num());
  std::stringstream ss;
  nested.serialize(ss);
  nested_[bin] = std::make_shared<RecursiveTable1D>(ss);
}

void RecursiveTable1D::insert(const std::vector<int>& bins,
    const RTable& nested) {
  ASSERT(static_cast<int>(bins.size()) == 1, "Error");
  insert(bins[0], nested);
}

double RecursiveTable1D::percent_nested() const {
  int num_nested = 0;
  for (const std::shared_ptr<RecursiveTable1D>& n : nested_) {
    if (n) {
      ++num_nested;
    }
  }
  return static_cast<double>(num_nested)/num();
}

double RecursiveTable1D::linear_interpolation(const double value0) const {
  int i0, i02;
  const double xd0 = table_xd_(value0, bin_spacing(), num(), &i0, &i02);
  TRACE("xd0 " << xd0);
  TRACE("i0 " << i0);
  TRACE("i02 " << i02);
  RecursiveTable1D * nested = nested_[i0].get();
  TRACE("nested? " << nested);
  if (nested) {
    TRACE("begin nested");
    return nested->linear_interpolation(xd0);
  }
  return c00_(xd0, i0, i02);
}

double RecursiveTable1D::linear_interpolation(
    const std::vector<double>& values) const {
  ASSERT(static_cast<int>(values.size()) == 1, "Error");
  return linear_interpolation(values[0]);
}

RecursiveTable1D RecursiveTable1D::combine(
    const std::vector<const RecursiveTable1D *>& rtables) const {
  const RecursiveTable1D * first = rtables[0];
  RecursiveTable1D rcombined({{"num", str(first->num())}});

  { DEBUG("first, combine the base tables");
    std::vector<const Table1D *> tables;
    for (const auto& rtable : rtables) {
      tables.push_back(rtable);
    }
    Table1D combined;
    combined = combined.combine(tables);
    rcombined.set_data(combined.data());
  }

  DEBUG("added all nested tables to combined table");
  std::stringstream ss;
  const int num = static_cast<int>(rtables.size());
  std::vector<const RecursiveTable1D *> ntabs(num);
  for (int i0 = 0; i0 < rcombined.num(); ++i0) {
    if (rtables[0]->nested_[i0]) {
      for (int itab = 0; itab < static_cast<int>(rtables.size()); ++itab) {
        ASSERT(rtables[itab]->nested_[i0], "err");
        ntabs[itab] = &(*rtables[itab]->nested_[i0]);
      }
      rcombined.nested_[i0] = std::make_shared<RecursiveTable1D>(rcombined.combine(ntabs));
    }
  }
  DEBUG("done");
  return rcombined;
}

double RecursiveTable1D::forward_difference_interpolation(const double value0) const {
  FATAL("not implemented.");
}

bool RecursiveTable1D::is_equal(const RecursiveTable1D& table) const {
  if (Table1D::is_equal(table)) {
    if (feasst::is_equal_fstobj(nested_, table.nested_)) {
      return true;
    }
  }
  return false;
}

RecursiveTable2D::RecursiveTable2D(argtype * args) : Table2D(args) {
  resize(num0(), num1(), &nested_);
}
RecursiveTable2D::RecursiveTable2D(argtype args) : RecursiveTable2D(&args) { feasst_check_all_used(args); }
RecursiveTable2D::~RecursiveTable2D() {}

void RecursiveTable2D::serialize(std::ostream& ostr) const {
  Table2D::serialize(ostr);
  feasst_serialize_version(6034, ostr);
  feasst_serialize(nested_, ostr);
}

RecursiveTable2D::RecursiveTable2D(std::istream& istr) : Table2D(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6034, "version: " << version);
//  feasst_deserialize(&nested_, istr);
  int dim0;
  istr >> dim0;
  nested_.resize(dim0);
  for (int index0 = 0; index0 < dim0; ++index0) {
    int dim1;
    istr >> dim1;
    auto * n1 = &nested_[index0];
    n1->resize(dim1);
    for (int index1 = 0; index1 < dim1; ++index1) {
      int existing;
      istr >> existing;
      if (existing != 0) {
        (*n1)[index1] = std::make_shared<RecursiveTable2D>(istr);
      }
    }
  }
}

void RecursiveTable2D::insert(const int bin0, const int bin1, const RTable& nested) {
  ASSERT(bin0 >= 0 && bin0 < num0(), "bin0:" << bin0 << " must be >0 and < num0:" << num0());
  ASSERT(bin1 >= 0 && bin1 < num1(), "bin1:" << bin1 << " must be >0 and < num1:" << num1());
  std::stringstream ss;
  nested.serialize(ss);
  nested_[bin0][bin1] = std::make_shared<RecursiveTable2D>(ss);
}

void RecursiveTable2D::insert(const std::vector<int>& bins,
    const RTable& nested) {
  ASSERT(static_cast<int>(bins.size()) == 2, "Error");
  insert(bins[0], bins[1], nested);
}

double RecursiveTable2D::linear_interpolation(const double value0, const double value1) const {
  int i0, i02, i1, i12;
  TRACE("value0 " << value0);
  TRACE("value1 " << value1);
  const double xd0 = table_xd_(value0, bin_spacing(0), num0(), &i0, &i02);
  TRACE("i0:" << i0 << " i02:" << i02);
  const double xd1 = table_xd_(value1, bin_spacing(1), num1(), &i1, &i12);
  TRACE("i1:" << i1 << " i12:" << i12);
  TRACE("xd0 " << xd0 << " xd1 " << xd1);
  RecursiveTable2D * nested = nested_[i0][i1].get();
  TRACE("nested? " << nested);
  if (nested) {
    TRACE("begin nested " << value0 << " " << value1);
    return nested->linear_interpolation(xd0, xd1);
  }
  return c00_(xd0, xd1, i0, i02, i1, i12);
}

int RecursiveTable2D::num_data() const {
  int num = num0()*num1();
  for (const std::vector<std::shared_ptr<RecursiveTable2D> >& ns : nested_) {
    for (const std::shared_ptr<RecursiveTable2D>& n : ns) {
      if (n) {
        num += n->num0()*n->num1();
      }
    }
  }
  return num;
}

double RecursiveTable2D::percent_nested() const {
  int num_nested = 0;
  for (const auto& ns : nested_) {
    for (const auto& n : ns) {
      if (n) {
        ++num_nested;
      }
    }
  }
  return static_cast<double>(num_nested)/num0()/num1();
}

double RecursiveTable2D::linear_interpolation(
    const std::vector<double>& values) const {
  ASSERT(static_cast<int>(values.size()) == 2, "Error");
  return linear_interpolation(values[0], values[1]);
}

RecursiveTable2D RecursiveTable2D::combine(
    const std::vector<const RecursiveTable2D *>& rtables) const {
  const RecursiveTable2D * first = rtables[0];
  RecursiveTable2D rcombined({
    {"num0", str(first->num0())},
    {"num1", str(first->num1())},
  });

  { DEBUG("first, combine the base tables");
    std::vector<const Table2D *> tables;
    for (const auto& rtable : rtables) {
      tables.push_back(rtable);
    }
    Table2D combined;
    combined = combined.combine(tables);
    rcombined.set_data(combined.data());
  }

  DEBUG("added all nested tables to combined table");
  std::stringstream ss;
  const int num = static_cast<int>(rtables.size());
  std::vector<const RecursiveTable2D *> ntabs(num);
  for (int i0 = 0; i0 < rcombined.num0(); ++i0) {
    for (int i1 = 0; i1 < rcombined.num1(); ++i1) {
      if (rtables[0]->nested_[i0][i1]) {
        for (int itab = 0; itab < static_cast<int>(rtables.size()); ++itab) {
          ASSERT(rtables[itab]->nested_[i0][i1], "err");
          ntabs[itab] = &(*rtables[itab]->nested_[i0][i1]);
        }
        rcombined.nested_[i0][i1] = std::make_shared<RecursiveTable2D>(rcombined.combine(ntabs));
      }
    }
  }
  DEBUG("done");
  return rcombined;
}

bool RecursiveTable2D::is_equal(const RecursiveTable2D& table) const {
  if (Table2D::is_equal(table)) {
    if (feasst::is_equal_fstobj(nested_, table.nested_)) {
      return true;
    }
  }
  return false;
}

RecursiveTable3D::RecursiveTable3D(argtype * args) : Table3D(args) {
  resize(num0(), num1(), num2(), &nested_);
}
RecursiveTable3D::RecursiveTable3D(argtype args) : RecursiveTable3D(&args) { feasst_check_all_used(args); }
RecursiveTable3D::~RecursiveTable3D() {}

void RecursiveTable3D::serialize(std::ostream& ostr) const {
  Table3D::serialize(ostr);
  feasst_serialize_version(6279, ostr);
  feasst_serialize(nested_, ostr);
}

RecursiveTable3D::RecursiveTable3D(std::istream& istr) : Table3D(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6279, "version: " << version);
  //feasst_deserialize(&nested_, istr);
  int dim0;
  istr >> dim0;
  nested_.resize(dim0);
  for (int index0 = 0; index0 < dim0; ++index0) {
    int dim1;
    istr >> dim1;
    auto * n1 = &nested_[index0];
    n1->resize(dim1);
    for (int index1 = 0; index1 < dim1; ++index1) {
      int dim2;
      istr >> dim2;
      auto * n2 = &(*n1)[index1];
      n2->resize(dim2);
      for (int index2 = 0; index2 < dim2; ++index2) {
        int existing;
        istr >> existing;
        if (existing != 0) {
          (*n2)[index2] = std::make_shared<RecursiveTable3D>(istr);
        }
      }
    }
  }
}

void RecursiveTable3D::insert(const int bin0, const int bin1, const int bin2,
    const RTable& nested) {
  ASSERT(bin0 >= 0 && bin0 < num0(), "bin0:" << bin0 << " must be >0 and < num0:" << num0());
  ASSERT(bin1 >= 0 && bin1 < num1(), "bin1:" << bin1 << " must be >0 and < num1:" << num1());
  ASSERT(bin2 >= 0 && bin2 < num2(), "bin2:" << bin2 << " must be >0 and < num2:" << num2());
  //TRACE("bin0 " << bin0 << " bin1 " << bin1 << " bin2 " << bin2);
  //TRACE("num0 " << num0() << " num1 " << num1() << " num2 " << num2());
  //WARN("For testing only.");
  std::stringstream ss;
  nested.serialize(ss);
  //TRACE("ss:" << ss.str());
  nested_[bin0][bin1][bin2] = std::make_shared<RecursiveTable3D>(ss);
  //nested_[bin0][bin1][bin2] = std::make_shared<RecursiveTable3D>(nested);
}

void RecursiveTable3D::insert(const std::vector<int>& bins,
    const RTable& nested) {
  ASSERT(static_cast<int>(bins.size()) == 3, "Error");
  insert(bins[0], bins[1], bins[2], nested);
}

double RecursiveTable3D::percent_nested() const {
  int num_nested = 0;
  for (const auto& ns2 : nested_) {
    for (const auto& ns : ns2) {
      for (const auto& n : ns) {
        if (n) {
          ++num_nested;
        }
      }
    }
  }
  return static_cast<double>(num_nested)/num0()/num1()/num2();
}

double RecursiveTable3D::linear_interpolation(const double value0, const double value1,
    const double value2) const {
  int i0, i02, i1, i12, i2, i22;
  const double xd0 = table_xd_(value0, bin_spacing(0), num0(), &i0, &i02);
  const double xd1 = table_xd_(value1, bin_spacing(1), num1(), &i1, &i12);
  const double xd2 = table_xd_(value2, bin_spacing(2), num2(), &i2, &i22);
  RecursiveTable3D * nested = nested_[i0][i1][i2].get();
  if (nested) {
    //TRACE("begin nested " << value0 << " " << value1);
    //TRACE("xd0 " << xd0 << " xd1 " << xd1);
    return nested->linear_interpolation(xd0, xd1, xd2);
  }
  return c00_(xd0, xd1, xd2, i0, i02, i1, i12, i2, i22);
}

double RecursiveTable3D::linear_interpolation(
    const std::vector<double>& values) const {
  ASSERT(static_cast<int>(values.size()) == 3, "Error");
  return linear_interpolation(values[0], values[1], values[2]);
}

RecursiveTable3D RecursiveTable3D::combine(
    const std::vector<const RecursiveTable3D *>& rtables) const {
  const RecursiveTable3D * first = rtables[0];
  RecursiveTable3D rcombined({
    {"num0", str(first->num0())},
    {"num1", str(first->num1())},
    {"num2", str(first->num2())},
  });

  { DEBUG("first, combine the base tables");
    std::vector<const Table3D *> tables;
    for (const auto& rtable : rtables) {
      tables.push_back(rtable);
    }
    Table3D combined;
    combined = combined.combine(tables);
    rcombined.set_data(combined.data());
  }

  DEBUG("added all nested tables to combined table");
  std::stringstream ss;
  const int num = static_cast<int>(rtables.size());
  std::vector<const RecursiveTable3D *> ntabs(num);
  for (int i0 = 0; i0 < rcombined.num0(); ++i0) {
    for (int i1 = 0; i1 < rcombined.num1(); ++i1) {
      for (int i2 = 0; i2 < rcombined.num2(); ++i2) {
        if (rtables[0]->nested_[i0][i1][i2]) {
          for (int itab = 0; itab < static_cast<int>(rtables.size()); ++itab) {
            ASSERT(rtables[itab]->nested_[i0][i1][i2], "err");
            ntabs[itab] = &(*rtables[itab]->nested_[i0][i1][i2]);
          }
          rcombined.nested_[i0][i1][i2] = std::make_shared<RecursiveTable3D>(rcombined.combine(ntabs));
        }
      }
    }
  }
  DEBUG("done");
  return rcombined;
}

bool RecursiveTable3D::is_equal(const RecursiveTable3D& table) const {
  if (Table3D::is_equal(table)) {
    if (feasst::is_equal_fstobj(nested_, table.nested_)) {
      return true;
    }
  }
  return false;
}

RecursiveTable5D::RecursiveTable5D(argtype * args) : Table5D(args) {
  resize(num0(), num1(), num2(), num3(), num4(), &nested_);
}
RecursiveTable5D::RecursiveTable5D(argtype args) : RecursiveTable5D(&args) { feasst_check_all_used(args); }
RecursiveTable5D::~RecursiveTable5D() {}

void RecursiveTable5D::serialize(std::ostream& ostr) const {
  Table5D::serialize(ostr);
  feasst_serialize_version(2346, ostr);
  feasst_serialize(nested_, ostr);
}

RecursiveTable5D::RecursiveTable5D(std::istream& istr) : Table5D(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2346, "version: " << version);
//  feasst_deserialize(&nested_, istr);
  int dim0;
  istr >> dim0;
  nested_.resize(dim0);
  for (int index0 = 0; index0 < dim0; ++index0) {
    int dim1;
    istr >> dim1;
    auto * n1 = &nested_[index0];
    n1->resize(dim1);
    for (int index1 = 0; index1 < dim1; ++index1) {
      int dim2;
      istr >> dim2;
      auto * n2 = &(*n1)[index1];
      n2->resize(dim2);
      for (int index2 = 0; index2 < dim2; ++index2) {
        int dim3;
        istr >> dim3;
        auto *n3 = &(*n2)[index2];
        n3->resize(dim3);
        for (int index3 = 0; index3 < dim3; ++index3) {
          int dim4;
          istr >> dim4;
          auto * n4 = &(*n3)[index3];
          n4->resize(dim4);
          for (int index4 = 0; index4 < dim4; ++index4) {
            int existing;
            istr >> existing;
            if (existing != 0) {
              (*n4)[index4] = std::make_shared<RecursiveTable5D>(istr);
            }
          }
        }
      }
    }
  }
}

void RecursiveTable5D::insert(const int bin0, const int bin1, const int bin2,
    const int bin3, const int bin4, const RTable& nested) {
  ASSERT(bin0 >= 0 && bin0 < num0(), "bin0:" << bin0 << " must be >0 and < num0:" << num0());
  ASSERT(bin1 >= 0 && bin1 < num1(), "bin1:" << bin1 << " must be >0 and < num1:" << num1());
  ASSERT(bin2 >= 0 && bin2 < num2(), "bin2:" << bin2 << " must be >0 and < num2:" << num2());
  ASSERT(bin3 >= 0 && bin3 < num3(), "bin3:" << bin3 << " must be >0 and < num3:" << num3());
  ASSERT(bin4 >= 0 && bin4 < num4(), "bin4:" << bin4 << " must be >0 and < num4:" << num4());
  std::stringstream ss;
  nested.serialize(ss);
  nested_[bin0][bin1][bin2][bin3][bin4] = std::make_shared<RecursiveTable5D>(ss);
}

void RecursiveTable5D::insert(const std::vector<int>& bins,
    const RTable& nested) {
  ASSERT(static_cast<int>(bins.size()) == 5, "bin size:" << bins.size());
  insert(bins[0], bins[1], bins[2], bins[3], bins[4], nested);
}

double RecursiveTable5D::percent_nested() const {
  int num_nested = 0;
  for (const auto& ns4 : nested_) {
    for (const auto& ns3 : ns4) {
      for (const auto& ns2 : ns3) {
        for (const auto& ns1 : ns2) {
          for (const auto& n : ns1) {
            if (n) {
              ++num_nested;
            }
          }
        }
      }
    }
  }
  return static_cast<double>(num_nested)/num0()/num1()/num2()/num3()/num4();
}

double RecursiveTable5D::linear_interpolation(const double value0, const double value1,
    const double value2, const double value3, const double value4) const {
  int i0, i02, i1, i12, i2, i22, i3, i32, i4, i42;
  const double xd0 = table_xd_(value0, bin_spacing(0), num0(), &i0, &i02);
  const double xd1 = table_xd_(value1, bin_spacing(1), num1(), &i1, &i12);
  const double xd2 = table_xd_(value2, bin_spacing(2), num2(), &i2, &i22);
  const double xd3 = table_xd_(value3, bin_spacing(3), num3(), &i3, &i32);
  const double xd4 = table_xd_(value4, bin_spacing(4), num4(), &i4, &i42);
  RecursiveTable5D * nested = nested_[i0][i1][i2][i3][i4].get();
  if (nested) {
    //TRACE("begin nested " << value0 << " " << value1);
    //TRACE("xd0 " << xd0 << " xd1 " << xd1);
    return nested->linear_interpolation(xd0, xd1, xd2, xd3, xd4);
  }
  return c00_(xd0, xd1, xd2, xd3, xd4, i0, i02, i1, i12, i2, i22, i3, i32, i4, i42);
}

RecursiveTable5D RecursiveTable5D::combine(
    const std::vector<const RecursiveTable5D *>& rtables) const {
  const RecursiveTable5D * first = rtables[0];
  RecursiveTable5D rcombined({
    {"num0", str(first->num0())},
    {"num1", str(first->num1())},
    {"num2", str(first->num2())},
    {"num3", str(first->num3())},
    {"num4", str(first->num4())},
  });

  { DEBUG("first, combine the base tables");
    std::vector<const Table5D *> tables;
    for (const auto& rtable : rtables) {
      tables.push_back(rtable);
    }
    Table5D combined;
    combined = combined.combine(tables);
    rcombined.set_data(combined.data());
  }

  DEBUG("added all nested tables to combined table");
  std::stringstream ss;
  const int num = static_cast<int>(rtables.size());
  std::vector<const RecursiveTable5D *> ntabs(num);
  for (int i0 = 0; i0 < rcombined.num0(); ++i0) {
    for (int i1 = 0; i1 < rcombined.num1(); ++i1) {
      for (int i2 = 0; i2 < rcombined.num2(); ++i2) {
        for (int i3 = 0; i3 < rcombined.num3(); ++i3) {
          for (int i4 = 0; i4 < rcombined.num4(); ++i4) {
            if (rtables[0]->nested_[i0][i1][i2][i3][i4]) {
              for (int itab = 0; itab < static_cast<int>(rtables.size()); ++itab) {
                ASSERT(rtables[itab]->nested_[i0][i1][i2][i3][i4], "err");
                ntabs[itab] = &(*rtables[itab]->nested_[i0][i1][i2][i3][i4]);
              }
              rcombined.nested_[i0][i1][i2][i3][i4] = std::make_shared<RecursiveTable5D>(rcombined.combine(ntabs));
            }
          }
        }
      }
    }
  }
  DEBUG("done");
  return rcombined;
}

bool RecursiveTable5D::is_equal(const RecursiveTable5D& table) const {
  if (Table5D::is_equal(table)) {
    if (feasst::is_equal_fstobj(nested_, table.nested_)) {
      return true;
    }
  }
  return false;
}

double RecursiveTable5D::linear_interpolation(
    const std::vector<double>& values) const {
  ASSERT(static_cast<int>(values.size()) == 5, "values: " << feasst_str(values));
  return linear_interpolation(values[0], values[1], values[2], values[3], values[4]);
}

const RecursiveTable5D& RecursiveTable5D::nested(const std::vector<int>& bins) {
  ASSERT(static_cast<int>(bins.size()) == 5, bins.size());
  return *nested_[bins[0]][bins[1]][bins[2]][bins[3]][bins[4]];
}

RecursiveTable6D::RecursiveTable6D(argtype * args) : Table6D(args) {
  resize(num0(), num1(), num2(), num3(), num4(), num5(), &nested_);
}
RecursiveTable6D::RecursiveTable6D(argtype args) : RecursiveTable6D(&args) { feasst_check_all_used(args); }
RecursiveTable6D::~RecursiveTable6D() {}

void RecursiveTable6D::serialize(std::ostream& ostr) const {
  Table6D::serialize(ostr);
  feasst_serialize_version(5478, ostr);
  feasst_serialize(nested_, ostr);
}

RecursiveTable6D::RecursiveTable6D(std::istream& istr) : Table6D(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5478, "version: " << version);
//  feasst_deserialize(&nested_, istr);
  int dim0;
  istr >> dim0;
  nested_.resize(dim0);
  for (int index0 = 0; index0 < dim0; ++index0) {
    int dim1;
    istr >> dim1;
    auto * n1 = &nested_[index0];
    n1->resize(dim1);
    for (int index1 = 0; index1 < dim1; ++index1) {
      int dim2;
      istr >> dim2;
      auto * n2 = &(*n1)[index1];
      n2->resize(dim2);
      for (int index2 = 0; index2 < dim2; ++index2) {
        int dim3;
        istr >> dim3;
        auto *n3 = &(*n2)[index2];
        n3->resize(dim3);
        for (int index3 = 0; index3 < dim3; ++index3) {
          int dim4;
          istr >> dim4;
          auto * n4 = &(*n3)[index3];
          n4->resize(dim4);
          for (int index4 = 0; index4 < dim4; ++index4) {
            int dim5;
            istr >> dim5;
            auto * n5 = &(*n4)[index4];
            n5->resize(dim5);
            for (int index5 = 0; index5 < dim5; ++index5) {
              int existing;
              istr >> existing;
              if (existing != 0) {
                (*n5)[index5] = std::make_shared<RecursiveTable6D>(istr);
              }
            }
          }
        }
      }
    }
  }
}

void RecursiveTable6D::insert(const int bin0, const int bin1, const int bin2,
    const int bin3, const int bin4, const int bin5, const RTable& nested) {
  ASSERT(bin0 >= 0 && bin0 < num0(), "bin0:" << bin0 << " must be >0 and < num0:" << num0());
  ASSERT(bin1 >= 0 && bin1 < num1(), "bin1:" << bin1 << " must be >0 and < num1:" << num1());
  ASSERT(bin2 >= 0 && bin2 < num2(), "bin2:" << bin2 << " must be >0 and < num2:" << num2());
  ASSERT(bin3 >= 0 && bin3 < num3(), "bin3:" << bin3 << " must be >0 and < num3:" << num3());
  ASSERT(bin4 >= 0 && bin4 < num4(), "bin4:" << bin4 << " must be >0 and < num4:" << num4());
  ASSERT(bin5 >= 0 && bin5 < num5(), "bin5:" << bin5 << " must be >0 and < num5:" << num5());
  std::stringstream ss;
  nested.serialize(ss);
  nested_[bin0][bin1][bin2][bin3][bin4][bin5] = std::make_shared<RecursiveTable6D>(ss);
}

void RecursiveTable6D::insert(const std::vector<int>& bins,
    const RTable& nested) {
  ASSERT(static_cast<int>(bins.size()) == 6, "Error");
  insert(bins[0], bins[1], bins[2], bins[3], bins[4], bins[5], nested);
}

double RecursiveTable6D::percent_nested() const {
  int num_nested = 0;
  for (const auto& ns5 : nested_) {
    for (const auto& ns4 : ns5) {
      for (const auto& ns3 : ns4) {
        for (const auto& ns2 : ns3) {
          for (const auto& ns1 : ns2) {
            for (const auto& n : ns1) {
              if (n) {
                ++num_nested;
              }
            }
          }
        }
      }
    }
  }
  return static_cast<double>(num_nested)/num0()/num1()/num2()/num3()/num4()/num5();
}

double RecursiveTable6D::linear_interpolation(const double value0, const double value1,
    const double value2, const double value3, const double value4,
    const double value5) const {
  int i0, i02, i1, i12, i2, i22, i3, i32, i4, i42, i5, i52;
  const double xd0 = table_xd_(value0, bin_spacing(0), num0(), &i0, &i02);
  const double xd1 = table_xd_(value1, bin_spacing(1), num1(), &i1, &i12);
  const double xd2 = table_xd_(value2, bin_spacing(2), num2(), &i2, &i22);
  const double xd3 = table_xd_(value3, bin_spacing(3), num3(), &i3, &i32);
  const double xd4 = table_xd_(value4, bin_spacing(4), num4(), &i4, &i42);
  const double xd5 = table_xd_(value5, bin_spacing(5), num5(), &i5, &i52);
  TRACE("xd0 " << xd0 << " xd1 " << xd1 << " xd2 " << xd2 << " xd3 " << xd3 <<
       " xd4 " << xd4 << " xd5 " << xd5);
  RecursiveTable6D * nested = nested_[i0][i1][i2][i3][i4][i5].get();
  if (nested) {
    TRACE("begin nested");
    return nested->linear_interpolation(xd0, xd1, xd2, xd3, xd4, xd5);
  }
  return c00_(xd0, xd1, xd2, xd3, xd4, xd5, i0, i02, i1, i12, i2, i22, i3, i32, i4, i42, i5, i52);
}

double RecursiveTable6D::linear_interpolation(
    const std::vector<double>& values) const {
  ASSERT(static_cast<int>(values.size()) == 6, "Error");
  return linear_interpolation(values[0], values[1], values[2], values[3], values[4], values[5]);
}

bool RecursiveTable6D::is_equal(const RecursiveTable6D& table) const {
  if (!Table6D::is_equal(table)) {
    INFO("base table unequal");
  } else {
    if (feasst::is_equal_fstobj(nested_, table.nested_)) {
      return true;
    } else {
      INFO("nested unequal");
    }
  }
  return false;
}

RecursiveTable6D RecursiveTable6D::combine(
    const std::vector<const RecursiveTable6D *>& rtables) const {
  const RecursiveTable6D * first = rtables[0];
  RecursiveTable6D rcombined({
    {"num0", str(first->num0())},
    {"num1", str(first->num1())},
    {"num2", str(first->num2())},
    {"num3", str(first->num3())},
    {"num4", str(first->num4())},
    {"num5", str(first->num5())},
  });

  { DEBUG("first, combine the base tables");
    std::vector<const Table6D *> tables;
    for (const auto& rtable : rtables) {
      tables.push_back(rtable);
    }
    Table6D combined;
    combined = combined.combine(tables);
    rcombined.set_data(combined.data());
  }

  DEBUG("added all nested tables to combined table");
  std::stringstream ss;
  const int num = static_cast<int>(rtables.size());
  std::vector<const RecursiveTable6D *> ntabs(num);
  for (int i0 = 0; i0 < rcombined.num0(); ++i0) {
    for (int i1 = 0; i1 < rcombined.num1(); ++i1) {
      for (int i2 = 0; i2 < rcombined.num2(); ++i2) {
        for (int i3 = 0; i3 < rcombined.num3(); ++i3) {
          for (int i4 = 0; i4 < rcombined.num4(); ++i4) {
            for (int i5 = 0; i5 < rcombined.num5(); ++i5) {
              if (rtables[0]->nested_[i0][i1][i2][i3][i4][i5]) {
                for (int itab = 0; itab < static_cast<int>(rtables.size()); ++itab) {
                  ASSERT(rtables[itab]->nested_[i0][i1][i2][i3][i4][i5], "err");
                  ntabs[itab] = &(*rtables[itab]->nested_[i0][i1][i2][i3][i4][i5]);
                }
                rcombined.nested_[i0][i1][i2][i3][i4][i5] = std::make_shared<RecursiveTable6D>(rcombined.combine(ntabs));
              }
            }
          }
        }
      }
    }
  }
  DEBUG("done");
  return rcombined;
}

}  // namespace feasst
