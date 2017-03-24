/**
 * \file
 *
 * \brief
 */

#include "./analyze_scatter_aniso.h"

/**
 * Constructor
 */
AnalyzeScatterAniso::AnalyzeScatterAniso(Space *space, Pair *pair)
  : Analyze(space, pair) {
}
AnalyzeScatterAniso::AnalyzeScatterAniso(Space *space,
  Pair *pair,
  const char* fileName)
  : Analyze(space, pair, fileName) {
}

