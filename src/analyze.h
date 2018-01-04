/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef ANALYZE_H_
#define ANALYZE_H_

#include "./pair.h"
#include "./criteria_wltmmc.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class WLTMMC;

/**
 * This is the base class used to write custom analysis code.
 * Operates on the Space and Pair class periodically and prints results.
 *
 * In order to create a custom Analysis, you can follow two similar procedures.
 *
 * First, you may define a custom Analysis code in the same file as
 * "int main()" for C++. See an example of this in the canonical LJ
 * <a href="testcase/1_lj_2_nvt-mc_README.html">test case</a>.
 *
 * Second, you may copy existing "analyze_*" files and replace the class name
 * and header guards (e.g. BASECLASS_DERIVED_H_).
 */
class Analyze : public BaseRandom {
 public:
  /// Constructor
  Analyze(Pair* pair,
    /**
     * allowed string key pairs (e.g., dictionary):
     *
     *  nFreq : Analyze every this many steps
     *
     *  - (default): see initFreq()
     *
     *  nFreqPrint : Print output every this many steps
     *
     *  - (default): see initFreqPrint()
     */
    const argtype &args = argtype());

  // HWH: Depreciated: Constructor
  Analyze(shared_ptr<Pair> pair, const argtype &args = argtype())
    : Analyze(pair.get(), args) {}

  /// This constructor is not often used, but its purpose is to initialize
  /// for MC interface before using reconstruct to set object pointers.
  Analyze() { defaultConstruction_(); }

  /// Initialize number of steps between each analysis.
  void initFreq(const int nfreq = 1) { nFreq_ = nfreq; }

  /// Return number of steps between each analysis.
  int nFreq() const { return nFreq_; }

  /// Initialize number of steps between each print to file.
  void initFreqPrint(const int nfreq = 1) { nFreqPrint_ = nfreq; }

  // Depreciated: HWH: old spelling
  void initPrintFreq(const int nfreq = 1) { nFreqPrint_ = nfreq; }

  /// Return number of steps between each print.
  int nFreqPrint() const { return nFreqPrint_; }

  /// Initialize production flag. 1 is on, 0 is off. Default is 1.
  virtual void initProduction(const int flag = 1) { production_ = flag; }

  /// Return production state.
  int production() const { return production_; }

  /// Initialize the name of the file to print results.
  void initFileName(const char* fileName) { fileName_.assign(fileName); }

  /// Append "chars" onto the file name where results are printed.
  void appendFileName(const char* chars) {
    if (!fileName_.empty()) fileName_.append(chars);
  }

  /// Perform the analysis and update the accumulators.
  virtual void update() { update(0); }

  /// Perform the analysis and update the accumulators.
  virtual void update(
    /// Specify the macrostate bin corresponding to CriteriaWLTMMC order param.
    const int iMacro
    ) { if (iMacro < -1) {} }

  /// Print the analysis to a file (default: restart file)
  virtual void write() {
    if (!fileName_.empty()) writeRestart(fileName_.c_str()); }

  /// Print the analysis to a file for each macrostate in CriteriaWLTMMC.
  virtual void write(CriteriaWLTMMC *c) {if (c == NULL) {} }

  /// Monkey patch to modify restart at run time for parallel restarting.
  //  NOTE to HWH: this is beyond scope of original intent of class
  virtual void modifyRestart(shared_ptr<WLTMMC> mc) { if (mc == NULL) {} }

  /// Write restart file.
  virtual void writeRestart(const char* fileName) {
    writeRestartBase(fileName);
  }
  void writeRestartBase(const char* fileName);

  /// Construct from restart file.
  Analyze(Pair* pair, const char* fileName);

  virtual ~Analyze() {}
  virtual Analyze* clone(Pair* pair) const;
  shared_ptr<Analyze> cloneShrPtr(Pair* pair) {
    return cloneImpl(pair); }

  /// Reset object pointers.
  void reconstruct(Pair *pair);

  /// Return pointer to space from pair.
  Space* space() { return pair_->space(); }

 protected:
  Pair *pair_;
  int nFreq_;        //!< frequency for analysis
  int nFreqPrint_;   //!< frequency for printing
  string fileName_;  //!< file name to print analysis
  int production_;   //!< set to 1 if in production

  // default constructor
  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Analyze> cloneImpl(Pair* pair) const;
};

/// Factory method
shared_ptr<Analyze> makeAnalyze(Pair* pair,
  const char* fileName);

/// Factory method
shared_ptr<Analyze> makeAnalyze(Pair* pair, const argtype &args = argtype());

/// Factory method
shared_ptr<Analyze> makeAnalyze();

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // ANALYZE_H_
