#include "./pair_round_square.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairRoundSquare::PairRoundSquare(Space* space,
  const double rCut)  //!< interaction cut-off distance
  : Pair(space, rCut) {
  className_.assign("PairRoundSquare");
}
PairRoundSquare::PairRoundSquare(Space* space,
  const char* fileName)
  : Pair(space, fileName) {
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

