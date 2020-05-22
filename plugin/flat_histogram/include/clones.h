
#ifndef FEASST_FLAT_HISTOGRAM_CLONES_H_
#define FEASST_FLAT_HISTOGRAM_CLONES_H_

#include <sstream>
#include <memory>
#include <vector>
#include "monte_carlo/include/monte_carlo.h"
#include "flat_histogram/include/window.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/flat_histogram.h"

namespace feasst {

/**
  Container for initializing, running and analyzing groups of FlatHistogram
  MonteCarlo simulations.
 */
class Clones {
 public:
  Clones() {}

  /// Add a MonteCarlo.
  void add(std::shared_ptr<MonteCarlo> mc) { clones_.push_back(mc); }

  // HWH this becomes too complicated with deep copies
  // HWH user functional creation of MonteCarlo is less complex
//  /// Create all clones via deep copy of given object and given windows.
//  void create(const MonteCarlo& mc, const Window& window);

  /// Return the number of clones.
  int num() const { return static_cast<int>(clones_.size()); }

  /// Return a read-only clone.
  const MonteCarlo& clone(const int index) const;

  /// Return a writable clone.
  MonteCarlo * get_clone(const int index);

  /// Return the clones.
  std::vector<std::shared_ptr<MonteCarlo> > get_clones() { return clones_; }

  /**
    Assuming the the first clone is already initialized, run the first clone
    until it reaches a macrostate the overlaps with the next clone, and exchange
    this configuration (e.g., bottom up initialization).
    Repeat until all clones have been initialized.

    args:
    - attempt_batch: perform this many attempts in a batch between checking
      for overlap (default: 1).
    - max_batch: maximum number of batches. Infinite if -1 (default: -1).
   */
  void initialize(const argtype& args = argtype());

  /// Same as above, but for only one overlapping pair of clones, given by
  /// upper_index and upper_index - 1.
  void initialize(const int upper_index,
    const argtype& args = argtype());

  /**
    Run until all clones are complete.
    If OMP is available, run the clones in parallel threads until all clones
    are complete.

    args:
    - omp_batch: If OMP, while waiting for other threads to complete, run this
      many trials in between checking for completion (default: 1).
   */
  void run_until_complete(const argtype& args = argtype());

  /**
    Combine the bottom up initialization and run until complete.
    With OMP, the first clone will run until finding overlap with the second.
    Once overlap is found, the first and second run in parallel while the second
    finds overlap with the third.
    This is repeated until all clones are running in parallel.
   */
  void initialize_and_run_until_complete(
    const argtype& init_args = argtype(),
    const argtype& run_args = argtype());

  /// Return the FlatHistogram of a given clone index.
  FlatHistogram flat_histogram(const int index) const;

  /// Return the LnProbability of all clones.
  LnProbability ln_prob(const argtype& args = argtype()) const;

  /// Serialize
  void serialize(std::ostream& ostr) const;

  /// Deserialize
  explicit Clones(std::istream& istr);

 private:
  std::vector<std::shared_ptr<MonteCarlo> > clones_;

  void run_until_complete_omp_(const argtype& run_args,
                               const bool init = false,
                               const argtype& init_args = argtype());
  void run_until_complete_serial_();
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_CLONES_H_
