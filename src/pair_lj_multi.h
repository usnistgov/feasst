/**
 * \file
 *
 * \brief multi-particle type implementation of lennard-jones pair-wise interaction
 *
 */
#ifndef PAIR_LJ_MULTI_H_
#define PAIR_LJ_MULTI_H_

#include "./pair_lj.h"

class PairLJMulti : public PairLJ {
 public:
  PairLJMulti(Space* space, const double rCut);
  PairLJMulti(Space* space, const char* fileName);
  ~PairLJMulti() {}
  virtual PairLJMulti* clone(Space* space) const {
    PairLJMulti* p = new PairLJMulti(*this); p->reconstruct(space); return p;
  }

  // defaults in constructor
  void defaultConstruction();

  /// write restart file
  virtual void writeRestart(const char* fileName);

  int initEnergy();     //!< function to calculate forces, given positions

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);

  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                 const int &jtype);

  /// potential energy and forces of all particles
  double allPartEnerForce(const int flag);

  /// inner loop for potential energy and forces of all particles
  void allPartEnerForceInner(const double &r2, const double &dx,
    const double &dy, const double &dz, const int &itype, const int &jtype,
    const int &iMol, const int &jMol);

  /// initialize cut and shifted potential
  void cutShift(const int flag);

  /// initialize linear force shift potential
  void linearShift(const int flag);

  /// set i-j shifting
  void cutShiftijset(const int itype, const int jtype, const int flag);
  void linearShiftijset(const int itype, const int jtype, const int flag);

  /// set i-j to wca
  void initWCA(const int itype, const int jtype);

  /// initialize long-range correcitons
  void initLRC();

  /// initialize exponential type
  // U_LJ(s=r/sig; alpha)/eps=4*(s^-2alpha - s^-alpha)
  //   type0 is 12-6: alpha=6
  //   type1 is 24-12: alpha=12
  //   type2 is 2alpha - alpha; alpha=16.6755
  //   type3 is 2alpha - alpha; alpha=50
  //   type4 is 2alpha - alpha; alpha=128
  //   type5 is 2alpha - alpha; alpha=24
  //   type6 is 2alpha - alpha; alpha=18
  void initExpType(const int type);

  /// initialize screened electrostatic interaction (Yukawa)
  ///  U(r) = A*exp(-kappa*r)/r
  void initScreenedElectro(const double A, const double K) {
    yukawa_ = 1; yukawaA_ = A; yukawaK_ = K;
  }
  void initScreenedElectro(const double A, const double K, const int yukawa) {
    yukawa_ = yukawa; yukawaA_ = A; yukawaK_ = K;
  }

  /// set the order parameter
  void setOrder(const double order);

  /// add a gaussian on the potential
  void addGaussian(const double height, const double position,
                   const double spread);

  /// set lambda parameter
  void setLambdaij(const double iType, const double jType, const double lambda);

  // return the lrc contribution of one particle
  double computeLRC(const int ipart);

  // read-only access to protected variables
  vector<double> rCutMax() const { return rCutMax_; }

 protected:
  /// potential energy shift by constant for i-j type interactions
  vector<vector<double> > peShiftij_;

  /// potential energy shift by linear term for i-j type interactions
  vector<vector<double> > peLinearShiftij_;

  /// precalculation for tail corrections (long range corrections)
  vector<vector<double> > lrcPreCalc_;

  //!< exponential type. type0 is 12-6. type1 is 24-12. 2: 33.351-16.6755
  //  3: 100-50. 4: 256-128
  int expType_;

  int yukawa_;          //!< turn on yukawa interactions if 1
  double yukawaA_;      //!< U(r) = A*exp(-K*r)/r
  double yukawaK_;      //!< U(r) = A*exp(-K*r)/r
  double alpha_;        //!< exponential parameter

  // lambda parameters
  vector<vector<double> > lambda_;       //!< DOI: 10.1021/ja802124e
  int lambdaFlag_;      //!< default: 0 (off)

  // gaussian parameters
  int gaussian_;        //!< flag for guassian interacitons
  vector<vector<double> > gausParam_;
};

#endif  // PAIR_LJ_MULTI_H_

