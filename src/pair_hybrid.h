/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef PAIR_HYBRID_H_
#define PAIR_HYBRID_H_

#include <vector>
#include <string>
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * This class calls multiple other pair classes.
 */
class PairHybrid : public Pair {
 public:
  /// Constructor
  /// @param rCut depreciated parameter with little meaning.
  PairHybrid(Space* space, const double rCut = 0);

  /// Add pair class.
  void addPair(Pair* pair) {pairVec_.push_back(pair); }

  /// Return a pair with index in order of the pairs added.
  Pair* getPair(const int iPair) const { return pairVec_[iPair]; }

  /// Select one pair so that all functions below operate on just that pair.
  void selectOne(const int iPair) {
    selected_.resize(1); selected_[0] = iPair; pairPrint_ = iPair;
  }

  void initEnergy();     //!< function to calculate forces, given positions

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);

  double peTot();   // total potential energy of system
  double vrTot();   // total virial of system

  /// stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);
  void update(const double de);

  /// add particle(s)
  void addPart();

  /// delete one particle
  void delPart(const int ipart);

  /// delete particles
  void delPart(const vector<int> mpart);

  /// print configuraiton
  int printxyz(const char* fileName, const int initFlag,
    const std::string comment = "") {
    return pairVec_[pairPrint_]->printxyz(fileName, initFlag, comment);
  }

  /// Identify a particle as non physical or non physical.
  void ipartNotPhysical(const int ipart);
  void ipartIsPhysical(const int ipart);
  void allPartPhysical();

  /// read only access to protected variables
  int nPairs() const { return static_cast<int>(pairVec_.size()); }

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  PairHybrid(Space* space, const char* fileName);

  ~PairHybrid();
  virtual PairHybrid* clone(Space* space) const;
  void reconstruct(Space* space);

 protected:
  vector<Pair*> pairVec_;   //!< vector of pointers to pairs
  int pairPrint_;       //!< index of pair to print
  int clone_;           //!< number of times the object has been cloned

  /// list of selected pairs to compure. If null, compute all
  vector<int> selected_;

  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairHybrid> makePairHybrid(Space* space, const double rCut = 0);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_HYBRID_H_
