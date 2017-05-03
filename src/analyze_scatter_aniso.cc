#include "./analyze_scatter_aniso.h"

namespace feasst {

AnalyzeScatterAniso::AnalyzeScatterAniso(Space *space, Pair *pair)
  : Analyze(space, pair) {
}
AnalyzeScatterAniso::AnalyzeScatterAniso(Space *space,
  Pair *pair,
  const char* fileName)
  : Analyze(space, pair, fileName) {
}

}  // namespace feasst

