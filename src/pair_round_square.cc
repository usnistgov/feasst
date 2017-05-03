#include "./pair_round_square.h"

namespace feasst {

PairRoundSquare::PairRoundSquare(Space* space,
  const double rCut)  //!< interaction cut-off distance
  : Pair(space, rCut) {
  className_.assign("PairRoundSquare");
}
PairRoundSquare::PairRoundSquare(Space* space,
  const char* fileName)
  : Pair(space, fileName) {
}

}  // namespace feasst

