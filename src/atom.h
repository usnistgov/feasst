/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef SRC_ATOM_H_
#define SRC_ATOM_H_

#include "base_random.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class Space;

/**
 * Define new properties per atom/site.
 */
class Atom : public BaseRandom {
 public:
  /// Constructor
  Atom();

  /// Initialize name identifier of property.
  void initName(const std::string name) { name_ = name; }

  /// Initialize data type of the property.
  void initDataType(const std::string type = "int");

  /// Set value
  void setVal(const int iatom, const int val) { intVal_[iatom] = val; }

  /// Return integer value
  int intVal(const int iatom) const { return intVal_[iatom]; }

  /// Return integer value
  void val(const int iatom, int * val) const { *val = intVal_[iatom]; }

  /// Return integer vector
  void vctr(vector<int> * vec) const { *vec = intVal_; }

  /// Return integer vector
  vector<int> intVal() const { return intVal_; }

  /// Return double precision value
  double dblVal(const int iatom) const { return dblVal_[iatom]; }

  /// Return integer value
  void val(const int iatom, double * val) const { *val = dblVal_[iatom]; }

  /// Return double vector
  void vctr(vector<double> * vec) const { *vec = dblVal_; }

  /// Return string value
  std::string strVal(const int iatom) const { return strVal_[iatom]; }

  /// Return string value
  void val(const int iatom, std::string * val) const { *val = strVal_[iatom]; }

  /// Return string vector
  void vctr(vector<std::string> * vec) const { *vec = strVal_; }

  /// Update due to addition of particle from space.
  virtual void addPart(const Space &space);

  /// Update due to deletion of particle from space.
  virtual void delPart(const int ipart);

  /// Update due to swap of site.
  void swap(const int previous, const int finl);

  /// Return the size of the data vector.
  unsigned int size() const;

  /// Write the restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  explicit Atom(const char* fileName);

  /// Return name of the group.
  std::string name() const { return name_; }

 protected:
  std::string type_;
  std::string name_;
  vector<int> intVal_;
  vector<double> dblVal_;
  vector<std::string> strVal_;

  /// Set value
  virtual void setVal_(const Space &space, const int iatom) {
    intVal_[iatom] = 0; }

  /// Set the default values during construction.
  void defaultConstruction_();
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // SRC_ATOM_H_
