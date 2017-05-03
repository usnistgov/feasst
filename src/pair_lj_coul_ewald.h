/**
 * \file
 *
 * \brief lennard-jones pairwise interactions with long range coulombic interactions treated by ewald
 *
 */
#ifndef PAIR_LJ_COUL_EWALD_H_
#define PAIR_LJ_COUL_EWALD_H_

#include "./pair.h"

namespace feasst {

class erftable {
 public:
  erftable();
  ~erftable() {}
  void init(const double alpha, const double rCut);
  double eval(const double x) const;
 private:
  vector<double> vtab_;
  int n_;
  double ds_;
};

class PairLJCoulEwald : public Pair {
 public:
  PairLJCoulEwald(Space* space, const double rCut);
  PairLJCoulEwald(Space* space, const char* fileName);
  ~PairLJCoulEwald();
  virtual PairLJCoulEwald* clone(Space* space) const {
    PairLJCoulEwald* p = new PairLJCoulEwald(*this);
    p->reconstruct(space);
    return p;
  }

  /// write restart file
  void writeRestart(const char* fileName);

  int initEnergy();     //!< function to calculate forces, given positions
  /// calculate reciprical (Fourier) space forces, given positions
  void forcesFrr();
  /// compute standard long range contributions of all particles
  void lrcConf();

  /// calculate interaction energy contribution of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);

  /// function to calculate real-space interaction energy contribution a subset
  //  of particles
  double multiPartEnerReal(const vector<int> mpart, const int flag);
  double multiPartEnerRealAtomCut(const vector<int> mpart, const int flag);

  /// compute fourier space contributions of multiple particles
  void multiPartEnerFrr(const vector<int> mpart, const int flag);

  /// stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change to the system
  void update(const vector<int> mpart, const int flag, const char* uptype);

  /// compute standard long range contributions of one particle ipart
  double lrcOne(const int ipart);

  void selfCorrect(vector<int> mpart);   //!< compute self interactions of mpart

  double peTot();   //!< total potential energy of system

  void delPart(const int ipart);   //!< delete one particle
  void delPart(const vector<int> mpart);   //!< delete particles
  void addPart();                       //!< add one particle

  double alpha;               //!< ewald damping parameter

  /// set maximum wave vector and compute self interactions
  void k2maxset(const int k2max);

  /// initialize bulk simulation of SPCE waters, atoms listed as oxygen then
  //  two hydrogens
  void initBulkSPCE();
  void initBulkSPCE(const double alphatmp, const int kmax);

  /// initialize kspace, number of wave vectors and screening parameters
  void initKSpace(const double alphatmp, const int k2max);

  /// initialize with LAMMPS data file
  void initLMPData(const string fileName);
  void initLMPData(const char* fileName) { initLMPData(string(fileName)); }

  /// initialize with LAMMPS data file
  void initJSONData(const string fileName) { ASSERT(fileName == "684358558679",
    "JSON not implemented for charges"); }

  /// check size of class variables
  void sizeCheck();

  /// read-only access of protected variables
  double peLJ() const { return peLJ_; }
  double peLJone() const { return peLJone_; }
  double peLRC() const { return peLRC_; }
  double peLRCone() const { return peLRCone_; }
  double peQReal() const { return peQReal_; }
  double peQRealone() const { return peQRealone_; }
  double peQFrr() const { return peQFrr_; }
  double peQFrrone() const { return peQFrrone_; }
  double peQFrrSelf() const { return peQFrrSelf_; }
  double peQFrrSelfone() const { return peQFrrSelfone_; }
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
  /// total potential energy from fourier space charge interactions
  double peQFrr_;
  /// total potential energy from fourier space charge interactions
  double deQFrr_;
  /// potential energy from one particle fourier space charge interactions
  double peQFrrone_;
  /// total potential energy from self interactions in fourier space
  double peQFrrSelf_;
  /// change in potential energy from self interactions in fourier space
  double deQFrrSelf_;
  /// potential energy from some self interactions in fourier space
  double peQFrrSelfone_;
  //!< fourier space wave vectors for real part of Ewald sum
  vector<double> eikrx_;
  vector<double> eikry_;
  vector<double> eikrz_;
  vector<double> eikix_;
  vector<double> eikiy_;
  vector<double> eikiz_;
  //!< new one particle fourier space wave vectors for real part of Ewald sum
  vector<double> eikrxnew_;
  vector<double> eikrynew_;
  vector<double> eikrznew_;
  //!< new particle fourier space wave vectors for imaginary part of Ewald sum
  vector<double> eikixnew_;
  vector<double> eikiynew_;
  vector<double> eikiznew_;

  vector<double> strucfacr_;               //!< structure factor
  vector<double> strucfacrnew_;               //!< structure factor
  vector<double> strucfacrold_;               //!< structure factor
  vector<double> strucfaci_;               //!< structure factor
  vector<double> strucfacinew_;               //!< structure factor
  vector<double> strucfaciold_;               //!< structure factor
  int kmax_;              //!< maximum wave vector magnitude cut off
  int k2max_;              //!< maximum wave vector squared magnitude cut off
  int kxmax_;              //!< maximum wave vector in particluar dimension
  int kymax_;              //!< maximum wave vector squared magnitude cut off
  int kzmax_;              //!< maximum wave vector squared magnitude cut off
  vector<double> kexp_;   //!< pre-computed wave-vector prefactor for Ewald sum
  vector<vector<double> > kvec_;   //!< pre-computed wave-vector for Ewald sum
  vector<int> k_;   //!< pre-computed k integers for Ewald sum

  erftable erft_;   //!< tabular error function
};

}  // namespace feasst

#endif  // PAIR_LJ_COUL_EWALD_H_

