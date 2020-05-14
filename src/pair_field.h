/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_FIELD_H_
#define PAIR_FIELD_H_

#include "./pair.h"

namespace feasst {

/**
 * Compute interactions of particles with a field.
 */
class PairField : virtual public Pair {
 public:
  /// Constructor.
  PairField(Space* space);

  /*
   * \return potential energy of a selection of particles multiPart, and store
   * this the Pair class variable peSRone
   */
  double multiPartEner(const vector < int > multiPart, const int flag = 0);

  /*
   * Compute the potential energy of all particles and store this in the Pair
   * class variable peTot
   */
  void initEnergy();

  // stores, restores or updates variables to avoid order recompute
  //   of entire configuration after every change
  virtual void update(const vector < int > mpart, const int flag,
                      const char* uptype);

  PairField(Space* space, const char* fileName);
  ~PairField() {}
  virtual PairField* clone(Space* space) const {
    PairField* p = new PairField(*this); p->reconstruct(space); return p;
  }

  virtual void writeRestart(const char* fileName);

 protected:
  double deWall_ = 0.;

  /**
   * Loop through a list of sites to compute the interaction of the particles
   * with a field.
   * A field interaction is based on the absolute position of the particles.
   */
  double fieldLoopSite_(
    /// the list of sites indicies must be sorted before (optimization)
    const vector<int> &siteList);

  /**
   * Compute the interaction between a site and a field.
   */
  virtual void fieldSite_(
    const int &siteType,  //!< type of site
    double * energy,      //!< energy of interaction
    double * force,       //!< force of interaction
    const double &x,      //!< x-dimension absolute position
    const double &y,      //!< y-dimension absolute position
    const double &z       //!< z-dimension absolute position
    ) { ASSERT(0, "not implemented"); }

  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairField> makePairField(std::shared_ptr<Space> space);

}  // namespace feasst

#endif  // PAIR_FIELD_H_
