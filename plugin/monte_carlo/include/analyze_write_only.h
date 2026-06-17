
#ifndef FEASST_MONTE_CARLO_ANALYZE_WRITE_ONLY_H_
#define FEASST_MONTE_CARLO_ANALYZE_WRITE_ONLY_H_

#include <vector>
#include <memory>
#include <string>
#include <map>
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  This Analyze does not perform updates.
  This class is for developers.
 */
class AnalyzeWriteOnly : public Analyze {
 public:
  explicit AnalyzeWriteOnly(argtype * args);

  void set_trials_per_update(const int trials) override;

  void set_trials_per(const int trials) { set_trials_per_write(trials); }

  explicit AnalyzeWriteOnly(std::istream& istr) : Analyze(istr) {}
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_ANALYZE_WRITE_ONLY_H_
