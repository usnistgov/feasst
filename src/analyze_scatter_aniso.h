/**
 * This file is a stub or placeholder for an experimental class that is not part of this release.
 */

#ifndef ANALYZE_SCATTER_ANISO_H_
#define ANALYZE_SCATTER_ANISO_H_

#include "./analyze.h"
#include "./shape.h"
#include "./pair_tabular.h"
#ifdef FFTW_
  #include <fftw3.h>
#endif  // FFTW_

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class AnalyzeScatterAniso : public Analyze {
 public:
  AnalyzeScatterAniso(Space *space, Pair *pair);
  AnalyzeScatterAniso(Space *space, Pair *pair, const char* fileName);
  ~AnalyzeScatterAniso() {}
  AnalyzeScatterAniso* clone(Space* space, Pair* pair) const {
    AnalyzeScatterAniso* a = new AnalyzeScatterAniso(*this);
    a->reconstruct(space, pair); return a;
  }
  shared_ptr<AnalyzeScatterAniso> cloneShrPtr(Space* space, Pair* pair) const {
    return(std::static_pointer_cast<AnalyzeScatterAniso, Analyze>(
      cloneImpl(space, pair)));
  }

  /// update analysis every nFreq
  void update() { update(0); }
  void update(const int iMacro) { if (iMacro < 1) {} }

 protected:
  // clone design pattern
  virtual shared_ptr<Analyze> cloneImpl(Space* space, Pair *pair) const {
    shared_ptr<AnalyzeScatterAniso> a = make_shared<AnalyzeScatterAniso>(*this);
    a->reconstruct(space, pair); return a;
  }
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // ANALYZE_SCATTER_ANISO_H_

