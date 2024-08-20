
#ifndef FEASST_STEPPERS_SEEK_ANALYZE_H_
#define FEASST_STEPPERS_SEEK_ANALYZE_H_

#include <string>
#include <vector>
#include "monte_carlo/include/analyze.h"

namespace feasst {

class MonteCarlo;

class AnalyzeData {
 public:
  virtual double get(const Analyze& analyze) const = 0;
  virtual ~AnalyzeData() {}
};

/// Obtain Accumulator average in Analyze.
class AccumulatorAverage : public AnalyzeData {
 public:
  double get(const Analyze& analyze) const override;
  virtual ~AccumulatorAverage() {}
};

/// Obtain Accumulator sum
class AccumulatorSum : public AnalyzeData {
 public:
  double get(const Analyze& analyze) const override;
  virtual ~AccumulatorSum() {}
};

/// Obtain Accumulator sum of squared
class AccumulatorSumOfSquared : public AnalyzeData {
 public:
  double get(const Analyze& analyze) const override;
  virtual ~AccumulatorSumOfSquared() {}
};

/// Obtain the average moment of Accumulator
class AccumulatorMoment : public AnalyzeData {
 public:
  AccumulatorMoment(const int moment = 0) { moment_ = moment; }
  double get(const Analyze& analyze) const override;
  virtual ~AccumulatorMoment() {}
 private:
  int moment_;
};

/**
  Find Analyze with class name.
 */
class SeekAnalyze {
 public:
  /**
    Return the indices, where the first is mc.analyze index.
    If inside AnalyzeFactory, the second is the index of the factory.
    Otherwise, the second index is -1.
    Only the first match is returned.
   */
  std::vector<int> index(const std::string class_name,
                         const MonteCarlo& mc) const;

  /// Return an Analyzer of given name.
  const Analyze& reference(const std::string class_name,
                           const MonteCarlo& mc) const;

  /**
    For multistate Analyze with given class_name,
    Return average Accumulator as a function of state.
   */
  std::vector<double> multistate_data(
    const std::string class_name,
    const MonteCarlo& mc,
    /// Optionally specify where to get data from Analyze
    const AnalyzeData& get = AccumulatorAverage()) const;
};

}  // namespace feasst

#endif  // FEASST_STEPPERS_SEEK_ANALYZE_H_
