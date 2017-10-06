#ifndef PAIR_LJ_COUL_EWALD_H_
#define PAIR_LJ_COUL_EWALD_H_

#include <string>
#include <vector>
#include "./pair.h"
#include "./table.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Lennard Jones pairwise interactions with long range coulombic interactions
 * treated by the Ewald summation.
 */
class PairLJCoulEwald : public Pair {
 public:
  /// Constructor
  PairLJCoulEwald(Space* space,
    /// interaction cut-off distance
    const double rCut);

  double alpha;               //!< ewald damping parameter

  /**
   * resize the wave vector arrays for a given maximum wave vector
   *  compute self interactions in fourier space
   */
  void k2maxset(
    /// maximum wave vector magnitude cut off
    const int k2max);

  /// initialize bulk simulation of SPCE waters, atoms listed as oxygen then
  //  two hydrogens
  void initBulkSPCE();

  /**
   * Initialize bulk simulation of SPCE waters, atoms listed as oxygen then two
   * hydrogens. Also set parameters.
   */
  void initBulkSPCE(const double alphatmp,  //!< Ewald alpha parameter
    /// Ewald maximum wave vector in a given dimension.
    const int kmax);

  /// initialize kspace, number of wave vectors and screening parameters
  void initKSpace(const double alphatmp, const int k2max);

  /// initialize with LAMMPS data file
  void initLMPData(const string fileName);
  void initLMPData(const char* fileName) { initLMPData(string(fileName)); }

  /// initialize with LAMMPS data file
  void initJSONData(const string fileName) { ASSERT(fileName == "684358558679",
    "JSON not implemented for charges"); }

  void initEnergy();     //!< function to calculate forces, given positions

  /**
   *  compute interaction contributions of multiple particles
   *  flag==0, energy of old configuration (no moves, ins or dels)
   *  flag==1, new configuration, same number of particles
   *  flag==2, old configuration, preparing to delete (same as 0)
   *  flag==3, just inserted particle
   */
  double multiPartEner(const vector<int> multiPart, const int flag);

  /// function to calculate real-space interaction energy contribution a subset
  //  of particles
  double multiPartEnerReal(const vector<int> mpart, const int flag);
  double multiPartEnerRealAtomCut(const vector<int> mpart, const int flag);

  /**
   * reciprical space energy of multiple particles, mpart.
   *  flag==0, energy of old configuration (no moves, ins or dels)
   *  flag==1, new configuration, same number of particles
   *  flag==2, old configuration, preparing to delete (same as 0)
   *  flag==3, just inserted particle
   */
  void multiPartEnerFrr(const vector<int> mpart, const int flag);

  /**
   * stores, restores or updates variables to avoid recompute of entire
   * configuration after every change
   *  flag==0, energy of old configuration (no moves, ins or dels)
   *  flag==1, new configuration, same number of particles
   *  flag==2, old configuration, preparing to delete (same as 0)
   *  flag==3, just inserted particle
   *  flag==5, computing entire configuration
   */
  void update(const vector<int> mpart,  //!< particles involved in move
    const int flag,   //!< type of move
    /// description of update type
    const char* uptype);

  double peTot();   //!< total potential energy of system

  void delPart(const int ipart);   // HWH depreciate: delete one particle
  void delPart(const vector<int> mpart);   //!< delete particles
  void addPart();                       //!< add particle(s)

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

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  PairLJCoulEwald(Space* space, const char* fileName);
  ~PairLJCoulEwald();
  virtual PairLJCoulEwald* clone(Space* space) const {
    PairLJCoulEwald* p = new PairLJCoulEwald(*this);
    p->reconstruct(space);
    return p;
  }

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

  /**
   * reciprical space forces, energy and virial
   *  fourier space sum over wave vectors for Ewald sum
   */
  void forcesFrr_();

  /// compute standard long range contributions of all particles
  void lrcConf_();

  /// compute standard long range contributions of one particle ipart
  double lrcOne_(const int ipart);

  /**
   * computes self interaction energies to compensate for Ewald Sum
   *  of mpart particles
   *  if mpart is null, updates self interaction energy for entire system
   */
  void selfCorrect_(vector<int> mpart);

  erftable erft_;   //!< tabular error function
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_LJ_COUL_EWALD_H_

