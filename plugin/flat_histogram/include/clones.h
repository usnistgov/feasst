
#ifndef FEASST_FLAT_HISTOGRAM_CLONES_H_
#define FEASST_FLAT_HISTOGRAM_CLONES_H_

#include <string>
#include <sstream>
#include <memory>
#include <vector>
#include "monte_carlo/include/monte_carlo.h"
#include "flat_histogram/include/window.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/flat_histogram.h"

namespace feasst {

//class Checkpoint;

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

//  /// Add a checkpoint.
//  void set(const std::shared_ptr<Checkpoint> checkpoint);

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
    - omp_batch: If OMP, for each this many steps, check for completion and
      write aggregate ln_prob (default: 1e6).
    - ln_prob_file: file name of aggregate ln_prob. If empty (default),
      do not write the file.
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

  /// Set the number of Criteria iterations of all clones.
  void set_num_iterations(const int iterations);

  /// Return the FlatHistogram of a given clone index.
  FlatHistogram flat_histogram(const int index) const;

  /// Return the LnProbability of all clones.
  LnProbability ln_prob(const argtype& args = argtype()) const;

  /// Serialize
  void serialize(std::ostream& ostr) const;

  /// Deserialize
  explicit Clones(std::istream& istr);

  std::string serialize() {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }
  Clones deserialize(const std::string str) {
    std::stringstream ss(str);
    return Clones(ss);
  }

 private:
  std::vector<std::shared_ptr<MonteCarlo> > clones_;
//  std::shared_ptr<Checkpoint> checkpoint_;

  void run_until_complete_omp_(const argtype& run_args,
                               const bool init = false,
                               const argtype& init_args = argtype());
  void run_until_complete_serial_();
};

/// Construct Clones from a vector of checkpoint file names.
std::shared_ptr<Clones> MakeClones(const std::vector<std::string> file_names);

/// Construct Clones indexed file names.
/// For example, if checkpoints of individual clones were named as follows:
/// checkpoint0.rst, checkpoint1.rst, checkpoint2.rst,
/// then use `auto clones = MakeClones("checkpoint, 3, 0, ".fst");`
std::shared_ptr<Clones> MakeClones(const std::string prepend,
  const int num,
  const int min = 0,
  const std::string append = "");

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_CLONES_H_
