#ifndef PAIR_SS_H_
#define PAIR_SS_H_

#include <vector>
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Soft sphere pair-wise interaction, \f$ U=eps(sig/r)^n \f$
 */
class PairSS : public Pair {
 public:
  /// Constructor
  PairSS(Space* space, const double rCut);

  /// Initialize the potential parameter exponent (default below).
  void initExponent(const int n = 12) { n_ = n; }

  /// Return the exponential parameter.
  double exponent() const { return n_; }

  /// Return the reduced second virial coefficient,
  /// b2~ = b2(beta eps)^(-3/n==12) as a reference
  double b2reduced();

  /// Write restart file.
  virtual void writeRestart(const char* fileName);

  /// Construct from restart file.
  PairSS(Space* space, const char* fileName);
  virtual ~PairSS() {}
  virtual PairSS* clone(Space* space) const {
    PairSS* p = new PairSS(*this); p->reconstruct(space); return p;
  }

  /// function to calculate forces, given positions
  virtual void initEnergy();

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                 const int &jtype);

  /// stores, restores or updates variables to avoid order recompute of entire
  /// configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);

 protected:
  double peSR_;
  double deSR_;
  double n_;    // potential parameter for power order

  // defaults in constructor
  void defaultConstruction_();
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_SS_H_
