/**
 * \file
 *
 * \brief lennard-jones pairwise interactions with long range coulombic interactions treated by ewald
 *
 */
#ifndef PAIR_LJ_COUL_H_
#define PAIR_LJ_COUL_H_

#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class PairLJCoul : public Pair {
 public:
  PairLJCoul(Space* space, const double rCut);
  PairLJCoul(Space* space, const char* fileName);
  ~PairLJCoul();
  virtual PairLJCoul* clone(Space* space) const {
    PairLJCoul* p = new PairLJCoul(*this);
    p->reconstruct(space);
    return p;
  }

  /// write restart file
  void writeRestart(const char* fileName);

  void initEnergy();     //!< function to calculate forces, given positions

  /// compute standard long range contributions of all particles
  void lrcConf();

  /// calculate interaction energy contribution of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);

  /// function to calculate real-space interaction energy contribution a subset
  //  of particles
  double multiPartEnerReal(const vector<int> mpart, const int flag);
  double multiPartEnerRealAtomCut(const vector<int> mpart, const int flag);

  /// stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change to the system
  void update(const vector<int> mpart, const int flag, const char* uptype);

  /// compute standard long range contributions of one particle ipart
  double lrcOne(const int ipart);

  double peTot();   //!< total potential energy of system

  /// initialize bulk simulation of SPCE waters, atoms listed as oxygen then
  //  two hydrogens
  void initBulkSPCE();
  void initBulkSPCE(const double alphatmp, const int kmax);

  /// initialize with LAMMPS data file
  void initLMPData(const string fileName);
  void initLMPData(const char* fileName) { initLMPData(string(fileName)); }

  /// initialize with LAMMPS data file
  void initJSONData(const string fileName) { ASSERT(fileName == "684358558679",
    "JSON not implemented for charges"); }

  /// check size of class variables
  void sizeCheck();

  /// hard sphere from COM
  double hardSphereCOM = 0.;

  /// read-only access of protected variables
  double peLJ() const { return peLJ_; }
  double peLJone() const { return peLJone_; }
  double peLRC() const { return peLRC_; }
  double peLRCone() const { return peLRCone_; }
  double peQReal() const { return peQReal_; }
  double peQRealone() const { return peQRealone_; }
  vector<double> q() const {return q_; }

 protected:
  vector<double> q_;      //!< particle charges
  double peLJ_;     //!< total potential energy from lennard-jones interactions
  double deLJ_;     //!< lennard jones potential energy change
  double peLJone_;  //!< lennard jones potential energy from subset of particles
  /// total potential energy from standard long range corrections
  double peLRC_;
  /// change in potential energy from standard long range corrections
  double deLRC_;
  /// potential energy from subset of particle standard long range corrections
  double peLRCone_;
  /// total potential energy from real space charge interactions
  double peQReal_;
  /// total potential energy from real space charge interactions
  double deQReal_;
  /// potential energy from one particle real space charge interactions
  double peQRealone_;
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_LJ_COUL_H_

