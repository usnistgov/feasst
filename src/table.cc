/**
 * \file
 *
 * \brief data table with interpolation
 */

#include "./table.h"

/**
 * Constructor
 */
Table::Table() {
  defaultConstruction();
}
Table::Table(const char* fileName) {
  defaultConstruction();
  tabDims_ = fstoi("tabDims", fileName);
  tabType_ = fstos("tabType", fileName);
  vector<int> tabsizes;
  vector<int> order;
  for (int i = 0; i < tabDims_; ++i) {
    stringstream ss;
    ss << "dimorder" << i;
    order.push_back(fstoi(ss.str().c_str(), fileName));
    ss.str("");
    ss << "dimn" << i;
    tabsizes.push_back(fstoi(ss.str().c_str(), fileName));
    vector<double> lims;
    ss.str("");
    ss << "dimmin" << i;
    lims.push_back(fstod(ss.str().c_str(), fileName));
    ss.str("");
    ss << "dimmax" << i;
    lims.push_back(fstod(ss.str().c_str(), fileName));
    tablim_.push_back(lims);
  }

  vector<int> orderInv(tabDims_);
  for (int i = 0; i < tabDims_; ++i) {
    orderInv[order[i]] = i;
  }

  // read table
  std::ifstream fs(fileName);
  std::string line;
  const int n = myProd(tabsizes);
  const int nLines = numLines(fileName);
  if (nLines > n) {
     for (int i = 0; i < nLines - n; ++i) getline(fs, line);
  }
  double val;
  if (tabDims_ == 1) {
    // precompute
    d0_ = (tablim_[0][1] - tablim_[0][0]) / static_cast<int>(tabsizes[0] - 1);
    tab1_.resize(tabsizes[0]);
    for (int dim0 = 0; dim0 < tabsizes[0]; ++dim0) {
      fs >> val;
      tab1_[dim0] = val;
      getline(fs, line);
    }
  } else if (tabDims_ == 3) {
    // precompute
    d0_ = (tablim_[0][1] - tablim_[0][0]) / static_cast<int>(tabsizes[0] - 1);
    d1_ = (tablim_[1][1] - tablim_[1][0]) / static_cast<int>(tabsizes[1] - 1);
    d2_ = (tablim_[2][1] - tablim_[2][0]) / static_cast<int>(tabsizes[2] - 1);
    vector<int> map(3);
    tab3_.resize(tabsizes[0], vector<vector<double> >(
                 tabsizes[1], vector<double>(tabsizes[2])));
    for (int dim2 = 0; dim2 < tabsizes[orderInv[2]]; ++dim2) {
      map[orderInv[2]] = dim2;
      for (int dim1 = 0; dim1 < tabsizes[orderInv[1]]; ++dim1) {
        map[orderInv[1]] = dim1;
        for (int dim0 = 0; dim0 < tabsizes[orderInv[0]]; ++dim0) {
          map[orderInv[0]] = dim0;
          fs >> val;
          tab3_[map[0]][map[1]][map[2]] = val;
          getline(fs, line);
        }
      }
    }
  } else if (tabDims_ == 4) {
    // precompute
    d0_ = (tablim_[0][1] - tablim_[0][0]) / static_cast<int>(tabsizes[0] - 1);
    d1_ = (tablim_[1][1] - tablim_[1][0]) / static_cast<int>(tabsizes[1] - 1);
    d2_ = (tablim_[2][1] - tablim_[2][0]) / static_cast<int>(tabsizes[2] - 1);
    d3_ = (tablim_[3][1] - tablim_[3][0]) / static_cast<int>(tabsizes[3] - 1);
    vector<int> map(4);
    tab4_.resize(tabsizes[0], vector<vector<vector<double> > >(
                 tabsizes[1], vector<vector<double> >(
                 tabsizes[2], vector<double>(tabsizes[3]))));
    for (int dim3 = 0; dim3 < tabsizes[orderInv[3]]; ++dim3) {
      map[orderInv[3]] = dim3;
      for (int dim2 = 0; dim2 < tabsizes[orderInv[2]]; ++dim2) {
        map[orderInv[2]] = dim2;
        for (int dim1 = 0; dim1 < tabsizes[orderInv[1]]; ++dim1) {
          map[orderInv[1]] = dim1;
          for (int dim0 = 0; dim0 < tabsizes[orderInv[0]]; ++dim0) {
            map[orderInv[0]] = dim0;
            fs >> val;
            tab4_[map[0]][map[1]][map[2]][map[3]] = val;
            getline(fs, line);
          }
        }
      }
    }
  } else if (tabDims_ == 5) {
    // precompute
    d0_ = (tablim_[0][1] - tablim_[0][0]) / static_cast<int>(tabsizes[0] - 1);
    d1_ = (tablim_[1][1] - tablim_[1][0]) / static_cast<int>(tabsizes[1] - 1);
    d2_ = (tablim_[2][1] - tablim_[2][0]) / static_cast<int>(tabsizes[2] - 1);
    d3_ = (tablim_[3][1] - tablim_[3][0]) / static_cast<int>(tabsizes[3] - 1);
    d4_ = (tablim_[4][1] - tablim_[4][0]) / static_cast<int>(tabsizes[4] - 1);
    vector<int> map(5);
    tab5_.resize(tabsizes[0], vector<vector<vector<vector<double> > > >(
      tabsizes[1], vector<vector<vector<double> > >(
      tabsizes[2], vector<vector<double> >(
      tabsizes[3], vector<double>(tabsizes[4])))));
    for (int dim4 = 0; dim4 < tabsizes[orderInv[4]]; ++dim4) {
      map[orderInv[4]] = dim4;
      for (int dim3 = 0; dim3 < tabsizes[orderInv[3]]; ++dim3) {
        map[orderInv[3]] = dim3;
        for (int dim2 = 0; dim2 < tabsizes[orderInv[2]]; ++dim2) {
          map[orderInv[2]] = dim2;
          for (int dim1 = 0; dim1 < tabsizes[orderInv[1]]; ++dim1) {
            map[orderInv[1]] = dim1;
            for (int dim0 = 0; dim0 < tabsizes[orderInv[0]]; ++dim0) {
              map[orderInv[0]] = dim0;
              fs >> val;
              tab5_[map[0]][map[1]][map[2]][map[3]][map[4]] = val;
              getline(fs, line);
            }
          }
        }
      }
    }
  } else if (tabDims_ == 6) {
    // precompute
    d0_ = (tablim_[0][1] - tablim_[0][0]) / static_cast<int>(tabsizes[0] - 1);
    d1_ = (tablim_[1][1] - tablim_[1][0]) / static_cast<int>(tabsizes[1] - 1);
    d2_ = (tablim_[2][1] - tablim_[2][0]) / static_cast<int>(tabsizes[2] - 1);
    d3_ = (tablim_[3][1] - tablim_[3][0]) / static_cast<int>(tabsizes[3] - 1);
    d4_ = (tablim_[4][1] - tablim_[4][0]) / static_cast<int>(tabsizes[4] - 1);
    d5_ = (tablim_[5][1] - tablim_[5][0]) / static_cast<int>(tabsizes[5] - 1);
    vector<int> map(6);
    tab6_.resize(
      tabsizes[0], vector<vector<vector<vector<vector<double> > > > >(
      tabsizes[1], vector<vector<vector<vector<double> > > >(
      tabsizes[2], vector<vector<vector<double> > >(
      tabsizes[3], vector<vector<double> >(
      tabsizes[4], vector<double>(tabsizes[5]))))));
    for (int dim5 = 0; dim5 < tabsizes[orderInv[5]]; ++dim5) {
      map[orderInv[5]] = dim5;
      for (int dim4 = 0; dim4 < tabsizes[orderInv[4]]; ++dim4) {
        map[orderInv[4]] = dim4;
        for (int dim3 = 0; dim3 < tabsizes[orderInv[3]]; ++dim3) {
          map[orderInv[3]] = dim3;
          for (int dim2 = 0; dim2 < tabsizes[orderInv[2]]; ++dim2) {
            map[orderInv[2]] = dim2;
            for (int dim1 = 0; dim1 < tabsizes[orderInv[1]]; ++dim1) {
              map[orderInv[1]] = dim1;
              for (int dim0 = 0; dim0 < tabsizes[orderInv[0]]; ++dim0) {
                map[orderInv[0]] = dim0;
                fs >> val;
                tab6_[map[0]][map[1]][map[2]][map[3]][map[4]][map[5]] = val;
                getline(fs, line);
              }
            }
          }
        }
      }
    }
  } else {
    cout << "error: Unrecognized tabdim(" << tabDims_ << ")" << endl;
    exit(0);
  }
  min_ = 0.;
}

/**
 * defaults in constructor
 */
void Table::defaultConstruction() {
}

/**
 * linear interpolation
 */
double Table::interpolate(const double val0) {
  const int i0 = (val0 - tablim_[0][0]) / d0_;
  int i02 = i0 + 1;
  if (i02 == static_cast<int>(tab1_.size())) i02 = i0;
  const double v0 = tablim_[0][0] + i0 * d0_, vv0 = v0 + d0_;
  const double xd0 = (val0 - v0) / (vv0 - v0);
  return tab1_[i0]*(1-xd0) + xd0*tab1_[i02];
}

/**
 * bilinear interpolation
 */
double Table::interpolate(const double val0, const double val1) {
  // cout << "d0 " << d0_ << " d1 " << d1_ << endl;
  const int i0 = (val0 - tablim_[0][0]) / d0_;
  const int i1 = (val1 - tablim_[1][0]) / d1_;
  // cout << "val0 " << val0 << " tablim00 " << tablim_[0][0] << endl;
  // cout << "val1 " << val1 << " tablim10 " << tablim_[1][0] << endl;
  // cout << "i0 " << i0 << " i1 " << i1 << endl;
  int i02 = i0 + 1, i12 = i1 + 1;
  if (i02 == static_cast<int>(tab2_.size())) i02 = i0;
  if (i12 == static_cast<int>(tab2_[0].size())) i12 = i1;
  // cout << "i02 " << i02 << " i12 " << i12 << endl;
  // cout << "static_cast<int>(tab2_[0].size()) "
  //   << static_cast<int>(tab2_[0].size()) << endl;
  const double v0 = tablim_[0][0] + i0 * d0_, vv0 = v0 + d0_;
  const double v1 = tablim_[1][0] + i1 * d1_, vv1 = v1 + d1_;
  // cout << "v0 " << v0 << " v1 " << v1 << endl;
  // cout << "vv0 " << vv0 << " vv1 " << vv1 << endl;
  const double xd0 = (val0 - v0) / (vv0 - v0);
  const double xd1 = (val1 - v1) / (vv1 - v1);
  // cout << "xd0 " << xd0 << " xd1 " << xd1 << endl;

  const double c0 = tab2_[i0][i1] *(1-xd0) + xd0*tab2_[i02][i1];
  const double c1 = tab2_[i0][i12]*(1-xd0) + xd0*tab2_[i02][i12];

  return c0*(1-xd1)+xd1*c1;
}

/**
 * trilinear interpolation
 */
double Table::interpolate(const double val0, const double val1,
  const double val2) {
  // cout << "d0 " << d0_ << " d1 " << d1_ << " d2 " << d2_ << endl;
  const int i0 = (val0 - tablim_[0][0]) / d0_;
  const int i1 = (val1 - tablim_[1][0]) / d1_;
  const int i2 = (val2 - tablim_[2][0]) / d2_;
  // cout << "val0 " << val0 << " tablim00 " << tablim_[0][0] << endl;
  // cout << "val1 " << val1 << " tablim10 " << tablim_[1][0] << endl;
  // cout << "val2 " << val2 << " tablim20 " << tablim_[2][0] << endl;
  // cout << "i0 " << i0 << " i1 " << i1 << " i2 " << i2 << endl;
  int i02 = i0 + 1, i12 = i1 + 1, i22 = i2 + 1;
  // if (i0 + 1 == tablim_[0][1]) i02 = 0;
  if (i02 == static_cast<int>(tab3_.size())) i02 = i0;
  if (i12 == static_cast<int>(tab3_[0].size())) i12 = i1;
  if (i22 == static_cast<int>(tab3_[0][0].size())) i22 = i2;
  // cout << "i02 " << i02 << " i12 " << i12 << " i22 " << i22 << endl;
  // cout << "tab3_.size() " << tab3_.size() << " tab3_[0].size() "
  // << tab3_[0].size() << " tab3_[0][0].size() " << tab3_[0][0].size() << endl;
  const double v0 = tablim_[0][0] + i0 * d0_, vv0 = v0 + d0_;
  const double v1 = tablim_[1][0] + i1 * d1_, vv1 = v1 + d1_;
  const double v2 = tablim_[2][0] + i2 * d2_, vv2 = v2 + d2_;
  // cout << "v0 " << v0 << " v1 " << v1 << " v2 " << v2 << endl;
  // cout << "vv0 " << vv0 << " vv1 " << vv1 << " vv2 " << vv2 << endl;
  const double xd0 = (val0 - v0) / (vv0 - v0);
  const double xd1 = (val1 - v1) / (vv1 - v1);
  const double xd2 = (val2 - v2) / (vv2 - v2);
  // cout << "xd0 " << xd0 << " xd1 " << xd1 << " xd2 " << xd2 << endl;

  const double c00 = tab3_[i0][i1][i2] *     (1-xd0) + xd0*tab3_[i02][i1][i2];
  const double c10 = tab3_[i0][i12][i2] *   (1-xd0) + xd0*tab3_[i02][i12][i2];
  const double c01 = tab3_[i0][i1][i22] *   (1-xd0) + xd0*tab3_[i02][i1][i22];
  const double c11 = tab3_[i0][i12][i22] * (1-xd0) + xd0*tab3_[i02][i12][i22];
  // cout << "c000 " << tab3_[i0][i1][i2] << " c010 " << tab3_[i0][i12][i2]
  //   << " c001 " << tab3_[i0][i1][i22] << " c011 " << tab3_[i0][i12][i22]
  //   << endl;
  // cout << "c100 " << tab3_[i02][i1][i2] << " c110 " << tab3_[i02][i12][i2]
  //   << " c101 " << tab3_[i02][i1][i22] << " c111 " << tab3_[i02][i12][i22]
  //   << endl;
  // cout << "c00 " << c00 << " c10 " << c10 << " c01 " << c01 << " c11 "
  //   << c11 << endl;

  const double c0 = c00 * (1-xd1) + xd1*c10;
  const double c1 = c01 * (1-xd1) + xd1*c11;
  // cout << "c0 " << c0 << " c1 " << c1 << endl;

  return c0*(1-xd2)+xd2*c1;
}

/**
 * 4D linear interpolation
 */
double Table::interpolate(const double val0, const double val1,
  const double val2, const double val3) {
  const int i0 = (val0 - tablim_[0][0]) / d0_;
  const int i1 = (val1 - tablim_[1][0]) / d1_;
  const int i2 = (val2 - tablim_[2][0]) / d2_;
  const int i3 = (val3 - tablim_[3][0]) / d3_;
  // cout << "val0 " << val0 << " tablim00 " << tablim_[0][0] << endl;
  // cout << "val1 " << val1 << " tablim10 " << tablim_[1][0] << endl;
  // cout << "val2 " << val2 << " tablim20 " << tablim_[2][0] << endl;
  // cout << "val3 " << val3 << " tablim30 " << tablim_[3][0] << endl;
  // cout << "i0 " << i0 << " i1 " << i1 << " i2 " << i2 << " i3 " << i3
  //    << endl;
  int i02 = i0 + 1, i12 = i1 + 1, i22 = i2 + 1, i32 = i3 + 1;
  // if (i0 + 1 == tablim_[0][1]) i02 = 0;
  if (i02 == static_cast<int>(tab4_.size())) i02 = i0;
  if (i12 == static_cast<int>(tab4_[0].size())) i12 = i1;
  if (i22 == static_cast<int>(tab4_[0][0].size())) i22 = i2;
  if (i32 == static_cast<int>(tab4_[0][0][0].size())) i32 = i3;
  // cout << "i02 " << i02 << " i12 " << i12 << " i22 " << i22 << endl;
  const double v0 = tablim_[0][0] + i0 * d0_, vv0 = v0 + d0_;
  const double v1 = tablim_[1][0] + i1 * d1_, vv1 = v1 + d1_;
  const double v2 = tablim_[2][0] + i2 * d2_, vv2 = v2 + d2_;
  const double v3 = tablim_[3][0] + i3 * d3_, vv3 = v3 + d3_;
  // cout << "v0 " << v0 << " v1 " << v1 << " v2 " << v2 << endl;
  // cout << "vv0 " << vv0 << " vv1 " << vv1 << " vv2 " << vv2 << endl;
  const double xd0 = (val0 - v0) / (vv0 - v0);
  const double xd1 = (val1 - v1) / (vv1 - v1);
  const double xd2 = (val2 - v2) / (vv2 - v2);
  const double xd3 = (val3 - v3) / (vv3 - v3);
  // cout << "xd0 " << xd0 << " xd1 " << xd1 << " xd2 " << xd2 << endl;

  const double c000 = tab4_[i0][i1 ][i2 ][i3 ] * (1-xd0)
               + xd0*tab4_[i02][i1 ][i2 ][i3 ];
  const double c100 = tab4_[i0][i12][i2 ][i3 ] * (1-xd0)
               + xd0*tab4_[i02][i12][i2 ][i3 ];
  const double c010 = tab4_[i0][i1 ][i22][i3 ] * (1-xd0)
               + xd0*tab4_[i02][i1 ][i22][i3 ];
  const double c110 = tab4_[i0][i12][i22][i3 ] * (1-xd0)
               + xd0*tab4_[i02][i12][i22][i3 ];
  const double c001 = tab4_[i0][i1 ][i2 ][i32] * (1-xd0)
               + xd0*tab4_[i02][i1 ][i2 ][i32];
  const double c101 = tab4_[i0][i12][i2 ][i32] * (1-xd0)
               + xd0*tab4_[i02][i12][i2 ][i32];
  const double c011 = tab4_[i0][i1 ][i22][i32] * (1-xd0)
               + xd0*tab4_[i02][i1 ][i22][i32];
  const double c111 = tab4_[i0][i12][i22][i32] * (1-xd0)
               + xd0*tab4_[i02][i12][i22][i32];

  const double c00 = c000 * (1-xd1) + xd1*c100;
  const double c10 = c010 * (1-xd1) + xd1*c110;
  const double c01 = c001 * (1-xd1) + xd1*c101;
  const double c11 = c011 * (1-xd1) + xd1*c111;

  const double c0 = c00 * (1-xd2) + xd2*c10;
  const double c1 = c01 * (1-xd2) + xd2*c11;

  return c0*(1-xd3)+xd3*c1;
}

/**
 * 5D linear interpolation
 */
double Table::interpolate(const double val0, const double val1,
  const double val2, const double val3, const double val4) {
  // cout << "d0 " << d0_ << " d1 " << d1_ << " d2 " << d2_ << " d3 " << d3_
  //   << " d4 " << d4_ << endl;
  const int i0 = (val0 - tablim_[0][0]) / d0_;
  const int i1 = (val1 - tablim_[1][0]) / d1_;
  const int i2 = (val2 - tablim_[2][0]) / d2_;
  const int i3 = (val3 - tablim_[3][0]) / d3_;
  const int i4 = (val4 - tablim_[4][0]) / d4_;
  // cout << "val0 " << val0 << " tablim00 " << tablim_[0][0] << " tablim01 "
  //   << tablim_[0][1] << endl;
  // cout << "val1 " << val1 << " tablim10 " << tablim_[1][0] << " tablim11 "
  //   << tablim_[1][1] << endl;
  // cout << "val2 " << val2 << " tablim20 " << tablim_[2][0] << " tablim21 "
  //   << tablim_[2][1] << endl;
  // cout << "val3 " << val3 << " tablim30 " << tablim_[3][0] << " tablim31 "
  //   << tablim_[3][1] << endl;
  // cout << "val4 " << val4 << " tablim40 " << tablim_[4][0] << " tablim41 "
  //   << tablim_[4][1] << endl;
  // cout << "i0 " << i0 << " i1 " << i1 << " i2 " << i2 << " i3 " << i3
  //   << " i4 " << i4 << endl;
  int i02 = i0 + 1, i12 = i1 + 1, i22 = i2 + 1, i32 = i3 + 1, i42 = i4 + 1;
  if (i02 == static_cast<int>(tab5_.size())) i02 = i0;
  if (i12 == static_cast<int>(tab5_[0].size())) i12 = i1;
  if (i22 == static_cast<int>(tab5_[0][0].size())) i22 = i2;
  if (i32 == static_cast<int>(tab5_[0][0][0].size())) i32 = i3;
  if (i42 == static_cast<int>(tab5_[0][0][0][0].size())) i42 = i4;
  // cout << "i02 " << i02 << " i12 " << i12 << " i22 " << i22 << " i32 "
  //   << i32 << " i42 " << i42 << endl;

  const double v0 = tablim_[0][0] + i0 * d0_, vv0 = v0 + d0_;
  const double v1 = tablim_[1][0] + i1 * d1_, vv1 = v1 + d1_;
  const double v2 = tablim_[2][0] + i2 * d2_, vv2 = v2 + d2_;
  const double v3 = tablim_[3][0] + i3 * d3_, vv3 = v3 + d3_;
  const double v4 = tablim_[4][0] + i4 * d4_, vv4 = v4 + d4_;

  const double xd0 = (val0 - v0) / (vv0 - v0);
  const double xd1 = (val1 - v1) / (vv1 - v1);
  const double xd2 = (val2 - v2) / (vv2 - v2);
  const double xd3 = (val3 - v3) / (vv3 - v3);
  const double xd4 = (val4 - v4) / (vv4 - v4);
  // cout << "xds " << xd0 << " " << xd1 << " " << xd2 << " " << xd3 << " "
  //   << xd4 << endl;

  const double c0000 = tab5_[i0][i1 ][i2 ][i3 ][i4 ]*(1-xd0)
                + xd0*tab5_[i02][i1 ][i2 ][i3 ][i4 ];
  const double c0001 = tab5_[i0][i1 ][i2 ][i3 ][i42]*(1-xd0)
                + xd0*tab5_[i02][i1 ][i2 ][i3 ][i42];
  const double c0010 = tab5_[i0][i1 ][i2 ][i32][i4 ]*(1-xd0)
                + xd0*tab5_[i02][i1 ][i2 ][i32][i4 ];
  const double c0011 = tab5_[i0][i1 ][i2 ][i32][i42]*(1-xd0)
                + xd0*tab5_[i02][i1 ][i2 ][i32][i42];
  const double c0100 = tab5_[i0][i1 ][i22][i3 ][i4 ]*(1-xd0)
                + xd0*tab5_[i02][i1 ][i22][i3 ][i4 ];
  const double c0101 = tab5_[i0][i1 ][i22][i3 ][i42]*(1-xd0)
                + xd0*tab5_[i02][i1 ][i22][i3 ][i42];
  const double c0110 = tab5_[i0][i1 ][i22][i32][i4 ]*(1-xd0)
                + xd0*tab5_[i02][i1 ][i22][i32][i4 ];
  const double c0111 = tab5_[i0][i1 ][i22][i32][i42]*(1-xd0)
                + xd0*tab5_[i02][i1 ][i22][i32][i42];
  const double c1000 = tab5_[i0][i12][i2 ][i3 ][i4 ]*(1-xd0)
                + xd0*tab5_[i02][i12][i2 ][i3 ][i4 ];
  const double c1001 = tab5_[i0][i12][i2 ][i3 ][i42]*(1-xd0)
                + xd0*tab5_[i02][i12][i2 ][i3 ][i42];
  const double c1010 = tab5_[i0][i12][i2 ][i32][i4 ]*(1-xd0)
                + xd0*tab5_[i02][i12][i2 ][i32][i4 ];
  const double c1011 = tab5_[i0][i12][i2 ][i32][i42]*(1-xd0)
                + xd0*tab5_[i02][i12][i2 ][i32][i42];
  const double c1100 = tab5_[i0][i12][i22][i3 ][i4 ]*(1-xd0)
                + xd0*tab5_[i02][i12][i22][i3 ][i4 ];
  const double c1101 = tab5_[i0][i12][i22][i3 ][i42]*(1-xd0)
                + xd0*tab5_[i02][i12][i22][i3 ][i42];
  const double c1110 = tab5_[i0][i12][i22][i32][i4 ]*(1-xd0)
                + xd0*tab5_[i02][i12][i22][i32][i4 ];
  const double c1111 = tab5_[i0][i12][i22][i32][i42]*(1-xd0)
                + xd0*tab5_[i02][i12][i22][i32][i42];
  // cout << " c0000 " << c0000 << endl;
  // cout << " c0001 " << c0001 << endl;

  const double c000 = c0000*(1-xd1) + xd1*c1000;
  const double c001 = c0001*(1-xd1) + xd1*c1001;
  const double c010 = c0010*(1-xd1) + xd1*c1010;
  const double c011 = c0011*(1-xd1) + xd1*c1011;
  const double c100 = c0100*(1-xd1) + xd1*c1100;
  const double c101 = c0101*(1-xd1) + xd1*c1101;
  const double c110 = c0110*(1-xd1) + xd1*c1110;
  const double c111 = c0111*(1-xd1) + xd1*c1111;

  const double c00 = c000 * (1-xd2) + xd2*c100;
  const double c10 = c010 * (1-xd2) + xd2*c110;
  const double c01 = c001 * (1-xd2) + xd2*c101;
  const double c11 = c011 * (1-xd2) + xd2*c111;

  const double c0 = c00 * (1-xd3) + xd3*c10;
  const double c1 = c01 * (1-xd3) + xd3*c11;

  return c0*(1-xd4)+xd4*c1;
}

/**
 * 6D linear interpolation
 */
double Table::interpolate(const double val0, const double val1,
  const double val2, const double val3, const double val4, const double val5) {
  // cout << "d0 " << d0_ << " d1 " << d1_ << " d2 " << d2_ << " d3 " << d3_
  //   << " d4 " << d4_ << " d5 " << d5_ << endl;
  const int i0 = (val0 - tablim_[0][0]) / d0_;
  const int i1 = (val1 - tablim_[1][0]) / d1_;
  const int i2 = (val2 - tablim_[2][0]) / d2_;
  const int i3 = (val3 - tablim_[3][0]) / d3_;
  const int i4 = (val4 - tablim_[4][0]) / d4_;
  const int i5 = (val5 - tablim_[5][0]) / d5_;
  // cout << "val0 " << val0 << " tablim00 " << tablim_[0][0] << endl;
  // cout << "val1 " << val1 << " tablim10 " << tablim_[1][0] << endl;
  // cout << "val2 " << val2 << " tablim20 " << tablim_[2][0] << endl;
  // cout << "val3 " << val3 << " tablim30 " << tablim_[3][0] << endl;
  // cout << "val4 " << val4 << " tablim40 " << tablim_[4][0] << endl;
  // cout << "val5 " << val5 << " tablim50 " << tablim_[5][0] << endl;
  // cout << "i0 " << i0 << " i1 " << i1 << " i2 " << i2 << " i3 " << i3
  //   << " i4 " << i4 << " i5 " << i5 << endl;
  int i02 = i0 + 1, i12 = i1 + 1, i22 = i2 + 1, i32 = i3 + 1, i42 = i4 + 1,
      i52 = i5 + 1;
  if (i02 == static_cast<int>(tab6_.size())) i02 = i0;
  if (i12 == static_cast<int>(tab6_[0].size())) i12 = i1;
  if (i22 == static_cast<int>(tab6_[0][0].size())) i22 = i2;
  if (i32 == static_cast<int>(tab6_[0][0][0].size())) i32 = i3;
  if (i42 == static_cast<int>(tab6_[0][0][0][0].size())) i42 = i4;
  if (i52 == static_cast<int>(tab6_[0][0][0][0][0].size())) i52 = i5;

  const double v0 = tablim_[0][0] + i0 * d0_, vv0 = v0 + d0_;
  const double v1 = tablim_[1][0] + i1 * d1_, vv1 = v1 + d1_;
  const double v2 = tablim_[2][0] + i2 * d2_, vv2 = v2 + d2_;
  const double v3 = tablim_[3][0] + i3 * d3_, vv3 = v3 + d3_;
  const double v4 = tablim_[4][0] + i4 * d4_, vv4 = v4 + d4_;
  const double v5 = tablim_[5][0] + i5 * d5_, vv5 = v5 + d5_;

  const double xd0 = (val0 - v0) / (vv0 - v0);
  const double xd1 = (val1 - v1) / (vv1 - v1);
  const double xd2 = (val2 - v2) / (vv2 - v2);
  const double xd3 = (val3 - v3) / (vv3 - v3);
  const double xd4 = (val4 - v4) / (vv4 - v4);
  const double xd5 = (val5 - v5) / (vv5 - v5);

  const double c00000 = tab6_[i0][i1 ][i2 ][i3 ][i4 ][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i2 ][i3 ][i4 ][i5 ];
  const double c00010 = tab6_[i0][i1 ][i2 ][i3 ][i42][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i2 ][i3 ][i42][i5 ];
  const double c00100 = tab6_[i0][i1 ][i2 ][i32][i4 ][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i2 ][i32][i4 ][i5 ];
  const double c00110 = tab6_[i0][i1 ][i2 ][i32][i42][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i2 ][i32][i42][i5 ];
  const double c01000 = tab6_[i0][i1 ][i22][i3 ][i4 ][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i22][i3 ][i4 ][i5 ];
  const double c01010 = tab6_[i0][i1 ][i22][i3 ][i42][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i22][i3 ][i42][i5 ];
  const double c01100 = tab6_[i0][i1 ][i22][i32][i4 ][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i22][i32][i4 ][i5 ];
  const double c01110 = tab6_[i0][i1 ][i22][i32][i42][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i22][i32][i42][i5 ];
  const double c10000 = tab6_[i0][i12][i2 ][i3 ][i4 ][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i12][i2 ][i3 ][i4 ][i5 ];
  const double c10010 = tab6_[i0][i12][i2 ][i3 ][i42][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i12][i2 ][i3 ][i42][i5 ];
  const double c10100 = tab6_[i0][i12][i2 ][i32][i4 ][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i12][i2 ][i32][i4 ][i5 ];
  const double c10110 = tab6_[i0][i12][i2 ][i32][i42][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i12][i2 ][i32][i42][i5 ];
  const double c11000 = tab6_[i0][i12][i22][i3 ][i4 ][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i12][i22][i3 ][i4 ][i5 ];
  const double c11010 = tab6_[i0][i12][i22][i3 ][i42][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i12][i22][i3 ][i42][i5 ];
  const double c11100 = tab6_[i0][i12][i22][i32][i4 ][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i12][i22][i32][i4 ][i5 ];
  const double c11110 = tab6_[i0][i12][i22][i32][i42][i5 ]*(1-xd0)
                 + xd0*tab6_[i02][i12][i22][i32][i42][i5 ];

  const double c00001 = tab6_[i0][i1 ][i2 ][i3 ][i4 ][i52]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i2 ][i3 ][i4 ][i52];
  const double c00011 = tab6_[i0][i1 ][i2 ][i3 ][i42][i52]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i2 ][i3 ][i42][i52];
  const double c00101 = tab6_[i0][i1 ][i2 ][i32][i4 ][i52]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i2 ][i32][i4 ][i52];
  const double c00111 = tab6_[i0][i1 ][i2 ][i32][i42][i52]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i2 ][i32][i42][i52];
  const double c01001 = tab6_[i0][i1 ][i22][i3 ][i4 ][i52]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i22][i3 ][i4 ][i52];
  const double c01011 = tab6_[i0][i1 ][i22][i3 ][i42][i52]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i22][i3 ][i42][i52];
  const double c01101 = tab6_[i0][i1 ][i22][i32][i4 ][i52]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i22][i32][i4 ][i52];
  const double c01111 = tab6_[i0][i1 ][i22][i32][i42][i52]*(1-xd0)
                 + xd0*tab6_[i02][i1 ][i22][i32][i42][i52];
  const double c10001 = tab6_[i0][i12][i2 ][i3 ][i4 ][i52]*(1-xd0)
                 + xd0*tab6_[i02][i12][i2 ][i3 ][i4 ][i52];
  const double c10011 = tab6_[i0][i12][i2 ][i3 ][i42][i52]*(1-xd0)
                 + xd0*tab6_[i02][i12][i2 ][i3 ][i42][i52];
  const double c10101 = tab6_[i0][i12][i2 ][i32][i4 ][i52]*(1-xd0)
                 + xd0*tab6_[i02][i12][i2 ][i32][i4 ][i52];
  const double c10111 = tab6_[i0][i12][i2 ][i32][i42][i52]*(1-xd0)
                 + xd0*tab6_[i02][i12][i2 ][i32][i42][i52];
  const double c11001 = tab6_[i0][i12][i22][i3 ][i4 ][i52]*(1-xd0)
                 + xd0*tab6_[i02][i12][i22][i3 ][i4 ][i52];
  const double c11011 = tab6_[i0][i12][i22][i3 ][i42][i52]*(1-xd0)
                 + xd0*tab6_[i02][i12][i22][i3 ][i42][i52];
  const double c11101 = tab6_[i0][i12][i22][i32][i4 ][i52]*(1-xd0)
                 + xd0*tab6_[i02][i12][i22][i32][i4 ][i52];
  const double c11111 = tab6_[i0][i12][i22][i32][i42][i52]*(1-xd0)
                 + xd0*tab6_[i02][i12][i22][i32][i42][i52];

  const double c0000 = c00000*(1-xd1) + xd1*c10000;
  const double c0010 = c00010*(1-xd1) + xd1*c10010;
  const double c0100 = c00100*(1-xd1) + xd1*c10100;
  const double c0110 = c00110*(1-xd1) + xd1*c10110;
  const double c1000 = c01000*(1-xd1) + xd1*c11000;
  const double c1010 = c01010*(1-xd1) + xd1*c11010;
  const double c1100 = c01100*(1-xd1) + xd1*c11100;
  const double c1110 = c01110*(1-xd1) + xd1*c11110;

  const double c0001 = c00001*(1-xd1) + xd1*c10001;
  const double c0011 = c00011*(1-xd1) + xd1*c10011;
  const double c0101 = c00101*(1-xd1) + xd1*c10101;
  const double c0111 = c00111*(1-xd1) + xd1*c10111;
  const double c1001 = c01001*(1-xd1) + xd1*c11001;
  const double c1011 = c01011*(1-xd1) + xd1*c11011;
  const double c1101 = c01101*(1-xd1) + xd1*c11101;
  const double c1111 = c01111*(1-xd1) + xd1*c11111;

  const double c000 = c0000 * (1-xd2) + xd2*c1000;
  const double c100 = c0100 * (1-xd2) + xd2*c1100;
  const double c010 = c0010 * (1-xd2) + xd2*c1010;
  const double c110 = c0110 * (1-xd2) + xd2*c1110;

  const double c001 = c0001 * (1-xd2) + xd2*c1001;
  const double c101 = c0101 * (1-xd2) + xd2*c1101;
  const double c011 = c0011 * (1-xd2) + xd2*c1011;
  const double c111 = c0111 * (1-xd2) + xd2*c1111;

  const double c00 = c000 * (1-xd3) + xd3*c100;
  const double c10 = c010 * (1-xd3) + xd3*c110;
  const double c01 = c001 * (1-xd3) + xd3*c101;
  const double c11 = c011 * (1-xd3) + xd3*c111;

  const double c0 = c00 * (1-xd4) + xd4*c10;
  const double c1 = c01 * (1-xd4) + xd4*c11;

  return c0*(1-xd5)+xd5*c1;
}

/**
 * compute minimum value in table
 */
double Table::compute_min() const {
  if (tabDims_ == 1) {
    return *std::min_element(tab1_.begin(), tab1_.begin()+tab1_.size());
  } else if (tabDims_ == 2) {
    vector<double> mins;
    for (unsigned int j = 0; j < tab2_.size(); ++j) {
      mins.push_back(*std::min_element(tab2_[j].begin(),
                     tab2_[j].begin() + tab2_[j].size()));
    }
    return *std::min_element(mins.begin(), mins.begin()+mins.size());
  } else if (tabDims_ == 3) {
    return myMinElement(tab3_);
  } else if (tabDims_ == 4) {
    return myMinElement(tab4_);
  } else if (tabDims_ == 5) {
    return myMinElement(tab5_);
  } else if (tabDims_ == 6) {
    return myMinElement(tab6_);
  } else {
    ASSERT(0, "error in table for compute_min, urecognized tabDims_("
      << tabDims_ << ")");
  }
  return -1;
}

/**
 * compute minimium values as a function of one table dimension, dim
 */
void Table::compute_min_compress1d(const int dim) {
  ASSERT(tabDims_ == 4, "ERROR: compute_min_compress1d in Table class cannot be"
    << "utilized with tabDims(" << tabDims_ << ").");
  ASSERT(dim == 0, "ERROR: compute_min_compress1d in Table class cannot be"
    << "utilized with dims(" << dim << ").");

  tab1_.resize(tab4_.size());
  for (unsigned int i = 0; i < tab4_.size(); ++i) {
    vector<double> mins2;
    for (unsigned int k = 0; k < tab4_[i].size(); ++k) {
      vector<double> mins;
      for (unsigned int j = 0; j < tab4_[i][k].size(); ++j) {
        mins.push_back(*std::min_element(tab4_[i][k][j].begin(),
                       tab4_[i][k][j].begin()+tab4_[i][k][j].size()));
      }
      mins2.push_back(*std::min_element(mins.begin(),
                      mins.begin()+mins.size()));
    }
    tab1_[i] = *std::min_element(mins2.begin(), mins2.begin()+mins2.size());
  }
}

/**
 * compute maximum value in table
 */
double Table::compute_max() const {
  if (tabDims_ == 3) {
    return myMaxElement(tab3_);
  } else if (tabDims_ == 4) {
    return myMaxElement(tab4_);
  } else if (tabDims_ == 5) {
    return myMaxElement(tab5_);
  } else {
    ASSERT(0, "error in table for compute_max, urecognized tabDims_("
      << tabDims_ << ")");
    exit(0);
  }
}

/**
 * print table in hdf5 format
 */
void Table::printHDF5(const char* fileName) {
#ifdef HDF5_
  const int tabDims_ = 3;
  const H5std_string  DATASET_NAME("Compressed_Data");
  const unsigned long long DIM0 = tab3_.size(),
    DIM1 = tab3_[0].size(),
    DIM2 = tab3_[0][0].size();

    hsize_t dims[tabDims_] = { DIM0, DIM1, DIM2 };  // dataset dimensions
    hsize_t chunk_dims[tabDims_];


    if (DIM0 > 20) {
      for (int dim = 0; dim < tabDims_; ++dim) {
        chunk_dims[dim] = 20;
      }
    } else {
      chunk_dims[0] = DIM0;
      chunk_dims[1] = DIM1;
      chunk_dims[2] = DIM2;
      // chunk_dims = { DIM0, DIM1, DIM2 };  // chunk dimensions
    }
    double *buf = new double[DIM0*DIM1*DIM2];

    // Try block to detect exceptions raised by any of the calls inside it
    try {

  // Turn off the auto-printing when failure occurs so that we can
  // handle the errors appropriately
  Exception::dontPrint();

  // Create a new file using the default property lists.
  H5File file(fileName, H5F_ACC_TRUNC);

  // Create the data space for the dataset.
  DataSpace *dataspace = new DataSpace(tabDims_, dims);

  // Modify dataset creation property to enable chunking
  DSetCreatPropList  *plist = new  DSetCreatPropList;
  plist->setChunk(tabDims_, chunk_dims);

  // Set ZLIB (DEFLATE) Compression using level 6.
  // To use SZIP compression comment out this line.
  plist->setDeflate(6);

  // Uncomment these lines to set SZIP Compression
  // unsigned szip_options_mask = H5_SZIP_NN_OPTION_MASK;
  // unsigned szip_pixels_per_block = 16;
  // plist->setSzip(szip_options_mask, szip_pixels_per_block);

  // Create the dataset.
  DataSet *dataset = new DataSet(file.createDataSet( DATASET_NAME,
                          PredType::NATIVE_DOUBLE, *dataspace, *plist) );

  for (unsigned int i = 0; i < DIM0; i++)
  for (unsigned int j = 0; j < DIM1; j++)
  for (unsigned int k = 0; k < DIM2; k++)
    buf[k+DIM2*(j+DIM1*i)] = tab3_[i][j][k];
//    vector[iz+nz*(iy+ny*ix)] = tab3_[ix][iy][iz];
//        buf[k+DIM2i][j][k] = tab3_[i][j][k];

  // Write data to dataset.
  dataset->write(buf, PredType::NATIVE_DOUBLE);

  // Close objects and file.  Either approach will close the HDF5 item.
  delete dataspace;
  delete dataset;
  delete plist;
  delete[] buf;
  file.close();

    }  // end of try block

    // catch failure caused by the H5File operations
    catch(FileIException error)
    {
  error.printError();
    }

    // catch failure caused by the DataSet operations
    catch(DataSetIException error)
    {
  error.printError();
    }

    // catch failure caused by the DataSpace operations
    catch(DataSpaceIException error)
    {
  error.printError();
    }
#endif  // HDF5_
}


