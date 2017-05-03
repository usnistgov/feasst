/**
 * \file
 *
 * \brief
 *
 */

#ifndef ANALYZEORIENT_H_
#define ANALYZEORIENT_H_

#include "./analyze.h"

namespace feasst {

class AnalyzeOrient : public Analyze {
 public:
  AnalyzeOrient(Space *space, Pair *pair);
  AnalyzeOrient(Space *space, Pair *pair, const char* fileName);
  ~AnalyzeOrient() {}
  AnalyzeOrient* clone(Space* space, Pair* pair) const {
    AnalyzeOrient* a = new AnalyzeOrient(*this);
    a->reconstruct(space, pair); return a;
  }
  shared_ptr<AnalyzeOrient> cloneShrPtr(Space* space, Pair* pair) const {
    return(std::static_pointer_cast<AnalyzeOrient, Analyze>(
      cloneImpl(space, pair)));
  }
  void defaultConstruction();
  void writeRestart(const char* fileName);

  /// update analysis every nFreq
  void update() { update(0); }
  void update(const int iMacro);

  /// print
  void print();
  void print(CriteriaWLTMMC *c) { print(); if (c == NULL) {} }

  double zbin;                   // histogram bin size

  /// read-only access
  vector<AccumulatorVec> zOrient() const { return zOrient_; }

 protected:
  Histogram h_;   //!< histogram to bin separation distances
  vector<AccumulatorVec> zOrient_;     // accumulator for average orientation

  // clone design pattern
  virtual shared_ptr<Analyze> cloneImpl(Space* space, Pair *pair) const {
    shared_ptr<AnalyzeOrient> a = make_shared<AnalyzeOrient>(*this)
    ; a->reconstruct(space, pair); return a;
  }
};

}  // namespace feasst

#endif  // ANALYZEORIENT_H_

