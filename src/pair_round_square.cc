/**
 * \file
 *
 * \brief 
 *
 */

#include "./pair_round_square.h"

/**
 * Constructor for pair_lj class requires the following
 */
PairRoundSquare::PairRoundSquare(Space* space,
  const double rCut)  //!< interaction cut-off distance
  : Pair(space, rCut) {
  className_.assign("PairRoundSquare");
}
PairRoundSquare::PairRoundSquare(Space* space,
  const char* fileName)
  : Pair(space, fileName) {
}

