/**
 * \file
 *
 * \brief ideal gas, no interactions
 *
 */
#ifndef PAIR_HYBRID_H_
#define PAIR_HYBRID_H_

#include "./pair.h"

namespace feasst {

class PairHybrid : public Pair {
 public:
  PairHybrid(Space* space, const double rCut);
  PairHybrid(Space* space, const char* fileName);
  ~PairHybrid();
  virtual PairHybrid* clone(Space* space) const;
  void reconstruct(Space* space);
  void writeRestart(const char* fileName);
  void defaultConstruction();

  int initEnergy();     //!< function to calculate forces, given positions

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);

  double peTot();   //!< total potential energy of system
  double vrTot();   //!< total virial of system

  /// stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);
  void update(const double de);

  /// add pair class
  void addPair(Pair* pair) {pairVec_.push_back(pair); }

  /// add particle(s)
  void addPart();

  /// delete one particle
  void delPart(const int ipart);

  /// delete particles
  void delPart(const vector<int> mpart);

  /// print configuraiton
  int printxyz(const char* fileName, const int initFlag,
    const std::string comment="") {
    return pairVec_[pairPrint_]->printxyz(fileName, initFlag, comment);
  }

  /// select one pair
  void selectOne(const int iPair) {
    selected_.resize(1); selected_[0] = iPair; pairPrint_ = iPair;
  }

  /// identify a particle as non physical or non physical
  void ipartNotPhysical(const int ipart);
  void ipartIsPhysical(const int ipart);
  void allPartPhysical();
  
  /// read only access to protected variables
  int nPairs() const { return static_cast<int>(pairVec_.size()); }

 protected:
  vector<Pair*> pairVec_;   //!< vector of pointers to pairs
  int pairPrint_;       //!< index of pair to print
  int clone_;           //!< number of times the object has been cloned

  /// list of selected pairs to compure. If null, compute all
  vector<int> selected_;
};

}  // namespace feasst

#endif  // PAIR_HYBRID_H_

