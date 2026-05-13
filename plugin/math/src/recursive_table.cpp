
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

void RecursiveTable1D::insert(const int bin, const RecursiveTable1D& nested) {
  ASSERT(bin >= 0 && bin < num(), "bin:" << bin << " must be >0 and < num:"
    << num());
  std::stringstream ss;
  nested.serialize(ss);
  nested_[bin] = std::make_shared<RecursiveTable1D>(ss);
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
  RecursiveTable1D * nested = nested_[i0].get();
  if (nested) {
    TRACE("begin nested");
    return nested->linear_interpolation(xd0);
  }
  return c00_(xd0, i0, i02);
}

double RecursiveTable1D::forward_difference_interpolation(const double value0) const {
  FATAL("not implemented.");
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

void RecursiveTable2D::insert(const int bin0, const int bin1, const RecursiveTable2D& nested) {
  ASSERT(bin0 >= 0 && bin0 < num0(), "bin0:" << bin0 << " must be >0 and < num0:" << num0());
  ASSERT(bin1 >= 0 && bin1 < num1(), "bin1:" << bin1 << " must be >0 and < num1:" << num1());
  std::stringstream ss;
  nested.serialize(ss);
  nested_[bin0][bin1] = std::make_shared<RecursiveTable2D>(ss);
}

double RecursiveTable2D::linear_interpolation(const double value0, const double value1) const {
  int i0, i02, i1, i12;
  const double xd0 = table_xd_(value0, bin_spacing(0), num0(), &i0, &i02);
  const double xd1 = table_xd_(value1, bin_spacing(1), num1(), &i1, &i12);
  RecursiveTable2D * nested = nested_[i0][i1].get();
  if (nested) {
    //INFO("begin nested " << value0 << " " << value1);
    //INFO("xd0 " << xd0 << " xd1 " << xd1);
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
  for (const std::vector<std::shared_ptr<RecursiveTable2D> >& ns : nested_) {
    for (const std::shared_ptr<RecursiveTable2D>& n : ns) {
      if (n) {
        ++num_nested;
      }
    }
  }
  return static_cast<double>(num_nested)/num0()/num1();
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
    const RecursiveTable3D& nested) {
  ASSERT(bin0 >= 0 && bin0 < num0(), "bin0:" << bin0 << " must be >0 and < num0:" << num0());
  ASSERT(bin1 >= 0 && bin1 < num1(), "bin1:" << bin1 << " must be >0 and < num1:" << num1());
  ASSERT(bin2 >= 0 && bin2 < num2(), "bin2:" << bin2 << " must be >0 and < num2:" << num2());
  //INFO("bin0 " << bin0 << " bin1 " << bin1 << " bin2 " << bin2);
  //INFO("num0 " << num0() << " num1 " << num1() << " num2 " << num2());
  //WARN("For testing only.");
  //std::stringstream ss2;
  //nested.serialize(ss2);
  //INFO("ss:" << ss2.str());
  //nested_[bin0][bin1][bin2] = std::make_shared<RecursiveTable3D>(ss2);
  nested_[bin0][bin1][bin2] = std::make_shared<RecursiveTable3D>(nested);
}

double RecursiveTable3D::percent_nested() const {
  int num_nested = 0;
  for (const std::vector<std::vector<std::shared_ptr<RecursiveTable3D> > >& ns2 : nested_) {
    for (const std::vector<std::shared_ptr<RecursiveTable3D> >& ns : ns2) {
      for (const std::shared_ptr<RecursiveTable3D>& n : ns) {
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
    //INFO("begin nested " << value0 << " " << value1);
    //INFO("xd0 " << xd0 << " xd1 " << xd1);
    return nested->linear_interpolation(xd0, xd1, xd2);
  }
  return c00_(xd0, xd1, xd2, i0, i02, i1, i12, i2, i22);
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
    const int bin3, const int bin4, const RecursiveTable5D& nested) {
  ASSERT(bin0 >= 0 && bin0 < num0(), "bin0:" << bin0 << " must be >0 and < num0:" << num0());
  ASSERT(bin1 >= 0 && bin1 < num1(), "bin1:" << bin1 << " must be >0 and < num1:" << num1());
  ASSERT(bin2 >= 0 && bin2 < num2(), "bin2:" << bin2 << " must be >0 and < num2:" << num2());
  ASSERT(bin3 >= 0 && bin3 < num3(), "bin3:" << bin3 << " must be >0 and < num3:" << num3());
  ASSERT(bin4 >= 0 && bin4 < num4(), "bin4:" << bin4 << " must be >0 and < num4:" << num4());
  std::stringstream ss;
  nested.serialize(ss);
  nested_[bin0][bin1][bin2][bin3][bin4] = std::make_shared<RecursiveTable5D>(ss);
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
    //INFO("begin nested " << value0 << " " << value1);
    //INFO("xd0 " << xd0 << " xd1 " << xd1);
    return nested->linear_interpolation(xd0, xd1, xd2, xd3, xd4);
  }
  return c00_(xd0, xd1, xd2, xd3, xd4, i0, i02, i1, i12, i2, i22, i3, i32, i4, i42);
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
    const int bin3, const int bin4, const int bin5, const RecursiveTable6D& nested) {
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

}  // namespace feasst
