#include "./analyze_scatter_aniso.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

AnalyzeScatterAniso::AnalyzeScatterAniso(Space *space, Pair *pair)
  : Analyze(space, pair) {
}
AnalyzeScatterAniso::AnalyzeScatterAniso(Space *space,
  Pair *pair,
  const char* fileName)
  : Analyze(space, pair, fileName) {
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

