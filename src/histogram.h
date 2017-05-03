/**
 * This class is used to compute histograms.
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <deque>
#include <iterator>
#include "./base.h"

namespace feasst {

class Histogram : public Base {
 public:
  //HWH: when is this empty constructor used? DEPRECIATE
  Histogram();  //!< Constructor
  explicit Histogram(const double binWidth);  //!< Constructor
  Histogram(const double binWidth, const int iType, const int jType);
  explicit Histogram(const char* fileName);
  virtual ~Histogram() {}
  Histogram* clone() const;
  shared_ptr<Histogram> cloneShrPtr() const;

  /// write restart file
  void writeRestart(const char* fileName);

  /// accumulate values
  virtual void accumulate(const double value);

  /// return bin
  double bin2m(const int bin) const { return min_ + (bin + 0.5)*binWidth_; }
  int bin(const double value) const {
    return myRound((value - bin2m(0))/binWidth_);
  }

  /// Center a bin on zero. Otherwise, zero lies on a boundary between bins.
  void centerZero();

  /// initialize the bin width
  void initBinWidth(double binWidth) { binWidth_ = binWidth; }

  /// print to file 
  void print(const char* fileName);

  /// count number of independent attempts to compute a histogram
  void count() { ++nNorm_; }

  /// return sum of elements in histogram
  double sum() const {
    return std::accumulate(histogram_.begin(), histogram_.end(), 0.);
  }

  /// return maximum element in histogram
  double maxElement() const {
    return *std::max_element(histogram_.begin(), histogram_.end());
  }

  /// return index of maximum element in histogram
  int maxElementBin() const;

  /// maximum value of the histogram
  double max() const { return max_; }
  
  /// minimum value of the histogram
  double min() const { return min_; }
  
  /// histogram data
  std::deque<double> hist() const { return histogram_; }
  
  /// size (or number of bins) of the histogram
  int size() const { return static_cast<int>(histogram_.size()); }
  
  /// width of bins in histogram is constant throughout the range
  double binWidth() const { return binWidth_; }
 
  // HWH: rename to nCount for connection with count()
  /// number of times histogram is computed. Used for normalization.
  long long nNorm() const { return nNorm_; }
  
  // HWH: when is this used? DEPRECIATE
  int iType() const { return iType_; }
  int jType() const { return jType_; }

 protected:
  double binWidth_;
  std::deque<double> histogram_;
  double max_;
  double min_;
  long long nNorm_;

  /// histograms may be described by pairs of types
  //  (e.g. radial distriubtion function)
  int iType_, jType_;

  /// flag to center histogram on zero (default, boundary is zero)
  int centerZero_;

  /// set the defaults in constructor
  void defaultConstruction_();

  // error messaging
  void mout_(const char* messageType, std::ostream& message) {
    myOut(messageType, message, className_, verbose_);
  }
};

}  // namespace feasst

#endif  // HISTOGRAM_H_

