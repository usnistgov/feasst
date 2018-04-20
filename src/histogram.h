/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <deque>
#include <iterator>
#include "./base.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * This class is used to compute one dimensional histograms with constant bin
 * width. By default, zero lies on a boundary between bins, but there is also
 * the option to center a histogram bin on zero.
 */
class Histogram : public Base {
 public:
  /**
   * Constructor requires binWidth, but if a non-physical value is employed,
   * then initBinWidth() must be utilized before accumulate().
   */
  explicit Histogram(const double binWidth = -1.);

  /// Constructor for checkpoint files.
  explicit Histogram(const char* fileName);

  virtual ~Histogram() {}

  /**
   * Clone raw pointer or copy of object which must be deleted to avoid a
   * memory leak.
   */
  Histogram* clone() const;

  /// Clone shared pointer.
  shared_ptr<Histogram> cloneShrPtr() const;

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Center a bin on zero. Otherwise, zero lies on a boundary between bins.
  void centerZero();

  /// Initialize the bin width.
  void initBinWidth(double binWidth) { binWidth_ = binWidth; }

  /// Increment the histogram by one for the bin corresponding to given value
  virtual void accumulate(const double value);

  /// Return value at center of bin, given bin.
  double bin2m(const int bin) const { return min_ + (bin + 0.5)*binWidth_; }

  /// Return bin, given value.
  int bin(const double value) const {
    return feasstRound((value - bin2m(0))/binWidth_);
  }

  /// Print to file.
  void write(const char* fileName);

  /// Count number of independent attempts to compute a histogram.
  void count() { ++nCount_; }

  /// Return number of times histogram is computed. Used for normalization.
  long long nCount() const { return nCount_; }

  /// Return sum of elements
  double sum() const {
    return std::accumulate(histogram_.begin(), histogram_.end(), 0.);
  }

  /// Return maximum element
  double maxElement() const {
    return *std::max_element(histogram_.begin(), histogram_.end());
  }

  /// Return index of maximum element
  int maxElementBin() const;

  /// Maximum value of the highest bin, e.g., the limit of the bins.
  double max() const { return max_; }

  /// Minimum value of the lowest bin, e.g., the limit of the bins.
  double min() const { return min_; }

  /// Return number of bins.
  int size() const { return static_cast<int>(histogram_.size()); }

  /// Return the width of the bins, which is constant throughout.
  double binWidth() const { return binWidth_; }

  /// Read-only access to raw histogram data.
  std::deque<double> hist() const { return histogram_; }

 protected:
  double binWidth_;
  std::deque<double> histogram_;
  double max_;
  double min_;
  long long nCount_;

  /// histograms may be described by pairs of types
  //  (e.g. radial distriubtion function)
  int iType_, jType_;

  /// flag to center histogram on zero (default, boundary is zero)
  int centerZero_;

  /// set the defaults in constructor
  void defaultConstruction_();
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // HISTOGRAM_H_

