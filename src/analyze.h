/**
 * \file
 *
 * Base class for analysis
 */

#ifndef ANALYZE_H_
#define ANALYZE_H_

#include "./space.h"
#include "./pair.h"
#include "./criteria_wltmmc.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class WLTMMC;

class Analyze : public BaseAll {
 public:
  Analyze() { defaultConstruction(); }
  Analyze(Space* space, Pair* pair);
  Analyze(Space *space, Pair* pair, const char* fileName);
  virtual ~Analyze() {}
  virtual Analyze* clone(Space* space, Pair* pair) const;

  // default constructor
  void defaultConstruction();

  // reset object pointers
  void reconstruct(Space* space, Pair *pair);

  /// write restart file
  virtual void writeRestart(const char* fileName) {
    writeRestartBase(fileName);
  }
  void writeRestartBase(const char* fileName);

  /// initialize frequencies and file names
  void initFreq(const int nfreq) { nFreq_ = nfreq; }
  void initPrintFreq(const int nfreq) { nFreqPrint_ = nfreq; }
  void initFileName(const char* fileName) { fileName_.assign(fileName); }

  /// append fileName
  void appendFileName(const char* chars) { fileName_.append(chars); }

  /// update analysis every nFreq
  virtual void update() { update(0); }
  virtual void update(const int iMacro) { if (iMacro < -1) {} }

  /// print analysis
  virtual void write() {
    if (!fileName_.empty()) writeRestart(fileName_.c_str()); }
  virtual void write(CriteriaWLTMMC *c) {if (c == NULL) {} }

  /// monkey patch to modify restart at run time
  //  NOTE to HWH: this is beyond scope of original intent of class
  virtual void modifyRestart(shared_ptr<WLTMMC> mc) { if (mc == NULL) {} }

  /// Initialize production.
  virtual void initProduction() { production_ = 1; }

  /// Initialize production flag. 1 is on, 0 is off. Default is 1.
  virtual void initProduction(const int flag) { production_ = flag; }

  // functions for read-only access of private data-members
  int nFreq() const { return nFreq_; }
  int nFreqPrint() const { return nFreqPrint_; }
  int production() const { return production_; }

 protected:
  Space *space_;
  Pair *pair_;
  int nFreq_;        //!< frequency for analysis
  int nFreqPrint_;   //!< frequency for printing
  string fileName_;  //!< file name to print analysis
  int production_;   //!< set to 1 if in production

  // clone design pattern
  virtual shared_ptr<Analyze> cloneImpl(Space* space, Pair* pair) const;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // ANALYZE_H_
