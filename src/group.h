/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef SRC_GROUP_H_
#define SRC_GROUP_H_

#include "atom.h"

namespace feasst {

class Space;

/**
 * Define a group of particles/sites for selection.
 */
class Group : public Atom {
 public:
  /// Constructor
  Group();

  /// Group by particle/molecule id.
  void molid(const int iMol);

//  /// Group by particle/molecule name.
//  void moltype(const char* name);
//  void moltype(const std::string name) { moltype(name.c_str()); }

  /// Return 1 if site is in group, 0 otherwise.
  int group(const int iAtom) const { return intVal_[iAtom]; }

  /// Return group vector
  vector<int> group() const { return intVal_; }

  /// Initialize Space which only contains particles in the group.
  void initSubSpace(shared_ptr<Space> space);

  /// Select a random molecule which is part of group.
  /// Molecule considered part of group if first atom is in the group.
  vector<int> randMol(Space * space);

  /// Write the restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  explicit Group(const char* fileName);

  // overload base class behavior
  void addPart(const Space &space);
  void delPart(const int ipart);

 private:
  /// Store criteria used for selection
  vector<int> molid_;
  vector<string> moltype_;

  /// Return 1 if iatom is part of group. Otherwise, 0.
  int iatomInGroup_(const Space &space, const int iatom);

  // Set value
  void setVal_(const Space &space, const int iatom) {
    intVal_[iatom] = iatomInGroup_(space, iatom);
  }

  // subspace
  shared_ptr<Space> subSpace_;
};

}  // namespace feasst

#endif  // SRC_GROUP_H_
