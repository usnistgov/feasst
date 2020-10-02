#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "steppers/include/seek_analyze.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

//int SeekAnalyze::index(const MonteCarlo& mc) const;
std::vector<int> SeekAnalyze::index(const std::string class_name,
                       const MonteCarlo& mc) const {
  DEBUG("Looking for " << class_name);
  DEBUG("num " << mc.num_analyzers());
  for (int index1 = 0; index1 < mc.num_analyzers(); ++index1) {
    const Analyze& an = mc.analyze(index1);
    std::string name = an.class_name();
    DEBUG("found " << name);
    if (name == class_name) {
      return {index1, -1};
    }
    if (name == "AnalyzeFactory") {
      const auto& ans = an.analyzers();
      for (int index2 = 0; index2 < static_cast<int>(ans.size()); ++index2) {
        name = ans[index2]->class_name();
        if (name == class_name) {
          return {index1, index2};
        }
      }
    }
  }
  return {-1, -1};
}

const Analyze& SeekAnalyze::reference(const std::string class_name,
    const MonteCarlo& mc) const {
  std::vector<int> ndx = index(class_name, mc);
  if (ndx[1] == -1) return mc.analyze(ndx[0]);
  if (ndx[0] != -1) return mc.analyze(ndx[0]).analyze(ndx[1]);
  FATAL(class_name << " not found.");
}

std::vector<double> SeekAnalyze::multistate_data(
    const std::string class_name,
    const MonteCarlo& mc,
    const AnalyzeData& get) const {
  const std::vector<int> ndx = index(class_name, mc);
  DEBUG(feasst_str(ndx));
  ASSERT(ndx[0] != -1, "class_name:" << class_name << " not found");
  ASSERT(ndx[1] != -1, "class_name:" << class_name << " is not multistate");
  std::vector<double> data;
  const auto& ans = mc.analyze(ndx[0]).analyzers();
  for (const auto& an : ans) data.push_back(get.get(*an));
  return data;
}

}  // namespace feasst
