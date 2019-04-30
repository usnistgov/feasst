
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_ACCUMULATOR_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_ACCUMULATOR_H_

#include <memory>
#include "core/include/accumulator.h"
#include "core/include/criteria.h"
#include "core/include/system.h"

namespace feasst {

/**
  HWH: Refactor this as an analyze class. Criteria base class should have
       state numbers.
  Accumulate a value for each bin of a macrostate.
 */
class MacrostateAccumulator {
 public:
  /// Accumulate a value.
  virtual void update(const System * system, const Criteria * criteria) = 0;

  /// Create a deep copy (of the derived class).
  std::shared_ptr<MacrostateAccumulator> deep_copy() { return deep_copy_impl_(); };

  /// Return the accumulator.
  const Accumulator& accumulator() const { return accumulator_; }

  /// Return the name for header.
  virtual std::string name() const = 0;

  virtual ~MacrostateAccumulator() {}

 protected:
  Accumulator accumulator_;

  virtual std::shared_ptr<MacrostateAccumulator> deep_copy_impl_() const = 0;
};

/**
  Apply a set of MacrostateAccumulator "types" to each bin of a macrostate.
 */
class MacrostateAccumulatorFactory {
 public:
  /// Add a new accumulator type for each bin.
  void add(std::shared_ptr<MacrostateAccumulator> tracker) { tracker_types_.push_back(tracker); }

  /// Resize the number of macrostates (or bins).
  void resize(const int size);

  /// Update the accumulators for a given bin.
  void update(const int bin, const System * system, const Criteria * criteria);

  /// Output the headers.
  std::string write_per_bin_header() const;

  /// Output the average and standard deviations.
  std::string write_per_bin(const int bin) const;

 private:
  std::vector<std::shared_ptr<MacrostateAccumulator> > tracker_types_;
  std::vector<std::vector<std::shared_ptr<MacrostateAccumulator> > > trackers_;
};

/**
  Store the system energy for each bin.
 */
class BinEnergy : public MacrostateAccumulator {
 public:
  void update(const System * system, const Criteria * criteria) override {
    accumulator_.accumulate(criteria->running_energy()); }
  std::string name() const override { return std::string("energy"); }
  virtual ~BinEnergy() {}

 protected:
  std::shared_ptr<MacrostateAccumulator> deep_copy_impl_() const override {
    return std::make_shared<BinEnergy>(*this); }
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_MACROSTATE_ACCUMULATOR_H_
