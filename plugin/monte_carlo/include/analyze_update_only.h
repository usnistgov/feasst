
#ifndef FEASST_MONTE_CARLO_ANALYZE_UPDATE_ONLY_H_
#define FEASST_MONTE_CARLO_ANALYZE_UPDATE_ONLY_H_

#include <vector>
#include <memory>
#include <string>
#include <map>
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  This Analyze does not perform writes.
  This class is for developers.
 */
class AnalyzeUpdateOnly : public Analyze {
 public:
  explicit AnalyzeUpdateOnly(argtype * args);

  void set_trials_per_write(const int trials) override;

  void set_trials_per(const int trials) { set_trials_per_update(trials); }

  explicit AnalyzeUpdateOnly(std::istream& istr) : Analyze(istr) {}
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_ANALYZE_UPDATE_ONLY_H_
