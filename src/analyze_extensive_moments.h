/**
 * \file analyze_extensive_moments.h
 *
 * \brief Store the extensive moments along the sampling order parameter of the system.  These are used to perform histogram extrapolation.
 *
 * This file is a stub or placeholder for an experimental class that is not part of this release.
 *
 */

#ifndef ANALYZEEXTMOMENTS_H_
#define ANALYZEEXTMOMENTS_H_

#include "./analyze.h"
#include "./criteria_wltmmc.h"
#include <iomanip>

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

#define EXTPRECISION 18 // Number of digits of precision to print extensive moments out as

class AnalyzeExtensiveMoments : public Analyze {
public:
    AnalyzeExtensiveMoments (Space *space, Pair *pair) : Analyze(space, pair) {}
    AnalyzeExtensiveMoments (Space *space, Pair *pair, const char* fileName) : Analyze(space, pair) {}
	~AnalyzeExtensiveMoments() {}
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif // ANALYZEEXTMOMENTS_H_
