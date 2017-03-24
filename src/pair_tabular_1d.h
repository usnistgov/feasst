/**
 * \file
 *
 * \brief tabular pair-wise interactions
 *
 *  This is only for istropic tabular potentials
 */
#ifndef PAIR_TABULAR_1D_H_
#define PAIR_TABULAR_1D_H_

#include "./pair.h"
#include "./table.h"

class PairTabular1D : public Pair {
 public:
  explicit PairTabular1D(Space* space);
  PairTabular1D(Space* space, const char* fileName);
  virtual ~PairTabular1D() {}
  virtual PairTabular1D* clone(Space* space) const {
    PairTabular1D* p = new PairTabular1D(*this);
    p->reconstruct(space);
    return p;
  }

  // defaults in constructor
  void defaultConstruction();

  /// write restart file
  virtual void writeRestart(const char* fileName);

  virtual int initEnergy();   //!< function to calculate forces, given positions

  /// potential energy of multiple particles
  double multiPartEner(const vector<int> multiPart, const int flag);

  /// potential energy and forces of all particles
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
    const int &jtype);

  /// inner loop for potential energy and forces of all particles
  double allPartEnerForce(const int flag);
  void allPartEnerForceInner(const double &r2, const double &dx,
    const double &dy, const double &dz, const int &itype, const int &jtype,
    const int &iMol, const int &jMol);

  /// stores, restores or updates variables to avoid order recompute of entire
  //  configuration after every change
  void update(const vector<int> mpart, const int flag, const char* uptype);

  /// delete one particle
  void delPart(const int ipart) {delPartBase(ipart); }
  void delPart(const vector<int> mpart) {delPartBase(mpart); }

  /// add one particle
  void addPart() {addPartBase(); }

  /// read tables
  void readTable(const char* fileName);

  /// read-only access of protected variables

 protected:
  vector<vector<double> > rCutInner_;      //!< hard sphere below inner cut-off
  double deSR_;             //!< lennard jones potential energy change
  vector<vector<shared_ptr<Table> > > tabij_;    //!< table for potential energy
  string tabFileName_;

  double tol_;      //!< table tolerance
};

#endif  // PAIR_TABULAR_1D_H_

