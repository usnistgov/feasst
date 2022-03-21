
#ifndef FEASST_FLAT_HISTOGRAM_COLLECT_MATRIX_SPLICE_H_
#define FEASST_FLAT_HISTOGRAM_COLLECT_MATRIX_SPLICE_H_

#include <string>
#include <sstream>
#include <memory>
#include <vector>
#include "monte_carlo/include/monte_carlo.h"
#include "flat_histogram/include/window.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/collection_matrix.h"
#include "steppers/include/seek_analyze.h"

namespace feasst {

class Checkpoint;
class Histogram;

/**
  Container for holding a group of FlatHistogram MonteCarlo simulations
  where each the collection matricies can be spliced together.
 */
class CollectionMatrixSplice {
 public:
  CollectionMatrixSplice() {}

  /// Add a MonteCarlo.
  void add(std::shared_ptr<MonteCarlo> mc) { clones_.push_back(mc); }

//  /// Create all clones via deep copy of given object and given windows.
//  CollectionMatrixSplice(const MonteCarlo& mc, const Window& window);

  /// Return the number of clones.
  int num() const { return static_cast<int>(clones_.size()); }

  /// Return a read-only clone.
  const MonteCarlo& clone(const int index) const;

  /// Return a writable clone.
  MonteCarlo * get_clone(const int index);

  /// Return the clones.
  std::vector<std::shared_ptr<MonteCarlo> > get_clones() { return clones_; }

  /// Return the FlatHistogram of a given clone index.
  FlatHistogram flat_histogram(const int index) const;

  /// Return the CollectionMatrix of a given clone index.
  TripleBandedCollectionMatrix collection_matrix(const int index) const;

  /// Loop through all clones and swap bounds based on number of interations.
  void adjust_bounds(const int min_window_size = 10);

  /// Return true if all clones are complete
  bool are_all_complete() const;

  /// Run for a number of hours.
  void run(const double hours);

  /// Return the complete collection matrix.
  TripleBandedCollectionMatrix collection_matrix() const;

  /// Return the complete probability distribution.
  LnProbability ln_prob() const;

  /// Serialize
  void serialize(std::ostream& ostr) const;

  /// Deserialize
  explicit CollectionMatrixSplice(std::istream& istr);

  std::string serialize() {
    std::stringstream ss;
    serialize(ss);
    return ss.str();
  }
  CollectionMatrixSplice deserialize(const std::string str) {
    std::stringstream ss(str);
    return CollectionMatrixSplice(ss);
  }

 private:
  std::vector<std::shared_ptr<MonteCarlo> > clones_;
  //std::shared_ptr<Checkpoint> checkpoint_;
};

/// Construct CollectionMatrixSplice
inline std::shared_ptr<CollectionMatrixSplice> MakeCollectionMatrixSplice() {
  return std::make_shared<CollectionMatrixSplice>(); }

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_COLLECT_MATRIX_SPLICE_H_
