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
  virtual void print() {}
  virtual void print(CriteriaWLTMMC *c) {if (c == NULL) {} }

  /// monkey patch to modify restart at run time
  //  NOTE to HWH: this is beyond scope of original intent of class
  virtual void modifyRestart(shared_ptr<WLTMMC> mc) { if (mc == NULL) {} }
  
  // functions for read-only access of private data-members
  int nFreq() const { return nFreq_; }
  int nFreqPrint() const { return nFreqPrint_; }

 protected:
  Space *space_;
  Pair *pair_;
  int nFreq_;        //!< frequency for analysis
  int nFreqPrint_;   //!< frequency for printing
  string fileName_;  //!< file name to print analysis

  // error messaging
  void mout_(const char* messageType, std::ostream& message) {
    myOut(messageType, message, className_, verbose_);
  }

  // clone design pattern
  virtual shared_ptr<Analyze> cloneImpl(Space* space, Pair* pair) const;
};

#endif  // ANALYZE_H_
