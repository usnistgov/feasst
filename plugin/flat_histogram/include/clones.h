
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
#include "steppers/include/seek_analyze.h"

namespace feasst {

class Checkpoint;
class Histogram;

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

  /// Add a checkpoint.
  void set(std::shared_ptr<Checkpoint> checkpoint);

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
  void initialize(argtype args = argtype());

  /// Same as above, but for only one overlapping pair of clones, given by
  /// upper_index and upper_index - 1.
  void initialize(const int upper_index, argtype args = argtype());

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
  void run_until_complete(argtype args = argtype());

  /**
    Combine the bottom up initialization and run until complete.
    With OMP, the first clone will run until finding overlap with the second.
    Once overlap is found, the first and second run in parallel while the second
    finds overlap with the third.
    This is repeated until all clones are running in parallel.
   */
  void initialize_and_run_until_complete(
    argtype run_args = argtype(),
    argtype init_args = argtype());

  /// Set the number of Criteria iterations of all clones.
  void set_num_iterations(const int iterations);

  /// Return the FlatHistogram of a given clone index.
  FlatHistogram flat_histogram(const int index) const;

  /// Stitch together and return the LnProbability of all clones.
  LnProbability ln_prob(
    /// Optionally return spliced macrostates, if not NULL.
    Histogram * macrostates = NULL,
    /// Optionally return spliced multistate data, it not NULL.
    std::vector<double> * multistate_data = NULL,
    /// Name of Analyze to extract data.
    const std::string analyze_name = "",
    /// Source of data in Analyze
    const AnalyzeData& get = AccumulatorAverage()) const;

  /// Same as above, but without ln_prob
  void stitch(
    Histogram * macrostates = NULL,
    std::vector<double> * multistate_data = NULL,
    const std::string analyze_name = "",
    const AnalyzeData& get = AccumulatorAverage()) const {
    ln_prob(macrostates, multistate_data, analyze_name, get); }

  /// Same as above, but without ln_prob or spliced macrostates.
  void stitch(std::vector<double> * multistate_data,
      const std::string analyze_name,
      const AnalyzeData& get) const {
    ln_prob(NULL, multistate_data, analyze_name, get); }

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
  std::shared_ptr<Checkpoint> checkpoint_;

  void run_until_complete_omp_(argtype run_args,
                               const bool init = false,
                               argtype init_args = argtype());
  void run_until_complete_serial_();
};

/// Construct Clones
inline std::shared_ptr<Clones> MakeClones() {
  return std::make_shared<Clones>(); }

/// Construct Clones from a vector of checkpoint file names.
std::shared_ptr<Clones> MakeClones(const std::vector<std::string> file_names);

/// Construct Clones indexed file names.
/// For example, if checkpoints of individual clones were named as follows:
/// checkpoint0.rst, checkpoint1.rst, checkpoint2.rst,
/// then use `auto clones = MakeClones("checkpoint, 3, 0, ".fst");`
std::shared_ptr<Clones> MakeClones(const std::string prepend,
  const int num,
  const int min = 0,
  const std::string append = ".fst");

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_CLONES_H_
