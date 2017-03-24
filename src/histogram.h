/**
 * \file
 *
 * \brief
 *
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <deque>
#include <iterator>
#include "./base_all.h"

class Histogram : public BaseAll {
 public:
  Histogram();
  explicit Histogram(const double binWidth);
  Histogram(const double binWidth, const int iType, const int jType);
  explicit Histogram(const char* fileName);
  virtual ~Histogram() {}
  Histogram* clone() const;
  shared_ptr<Histogram> cloneShrPtr() const;

  // defaults in constructor
  void defaultConstruction();

  /// write restart file
  void writeRestart(const char* fileName);

  /// accumulate values
  virtual void accumulate(const double value);

  /// return bin
  double bin2m(const int bin) const { return min_ + (bin + 0.5)*binWidth_; }
  int bin(const double value) const {
    return myRound((value - bin2m(0))/binWidth_);
  }

  /// count number of independent attempts to compute a histogram for
  //  normalization
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

  /// center the histogram on zero
  void centerZero();

  /// initialize the bin width
  void initBinWidth(double binWidth) { binWidth_ = binWidth; }

  /// print to file 
  void print(const char* fileName);

  /// read-only access of private data-members
  double max() const { return max_; }
  double min() const { return min_; }
  std::deque<double> hist() const { return histogram_; }
  int size() const { return static_cast<int>(histogram_.size()); }
  double binWidth() const { return binWidth_; }
  long long nNorm() const { return nNorm_; }
  int iType() const { return iType_; }
  int jType() const { return jType_; }

 protected:
  double binWidth_;           //!< bin width
  std::deque<double> histogram_;    //!< double-ended queue
  double max_;                      //!< maximum value in deque
  double min_;                      //!< minimum value in deque

  /// count number of times histogram is computed for normalization
  long long nNorm_;

  /// histograms may be described by pairs of types
  //  (e.g. radial distriubtion function)
  int iType_, jType_;

  /// flag to center histogram on zero (default, boundary is zero)
  int centerZero_;

  // error messaging
  void mout_(const char* messageType, std::ostream& message) {
    myOut(messageType, message, className_, verbose_);
  }
};

#endif  // HISTOGRAM_H_

