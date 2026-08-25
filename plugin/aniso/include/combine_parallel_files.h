
#ifndef FEASST_ANISO_COMBINE_PARALLEL_FILES_H_
#define FEASST_ANISO_COMBINE_PARALLEL_FILES_H_

#include <memory>
#include <string>
#include <vector>
#include <cstdint>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  This class is in active development and not recommended for general users.
  This class is used to parallelize BuildRecursiveTable.
 */
class CombineParallelFiles : public Action {
 public:
  //@{
  /** @name Arguments
    - type: type of file, which can be recursive_table or trained.
    - input_prefix: The beginning of the input file names for the multiple
      files to combine into one, before the processor index.
      The combined file is output as the input prefix and suffix without
      processor index.
    - input_suffix: The end of the input file after the processor index
      (default: empty).
    - num_processors: number of processors.
    - cleanup: If true (default), remove input files after combining.
   */
  explicit CombineParallelFiles(argtype args = argtype());
  explicit CombineParallelFiles(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override {}
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<CombineParallelFiles>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<CombineParallelFiles>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit CombineParallelFiles(std::istream& istr);
  virtual ~CombineParallelFiles();
  //@}
};

}  // namespace feasst

#endif  // FEASST_ANISO_COMBINE_PARALLEL_FILES_H_
