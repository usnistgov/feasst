/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./atom.h"
#include "./space.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Atom::Atom() {
  defaultConstruction_();
}

void Atom::defaultConstruction_() {
  initDataType();
  name_ = randomHash(5);
}

void Atom::addPart(const Space &space) {
  const int nAtomPrev = intVal_.size();
  if (type_ == "int") {
    for (int iPart = nAtomPrev; iPart < space.natom(); ++iPart) {
      intVal_.push_back(0);
      setVal_(space, iPart);
    }
  }
}

void Atom::delPart(const int ipart) {
  ASSERT(ipart < static_cast<int>(intVal_.size()),
    "cannot delete ipart(" << ipart << ") when "
    << "per atom quantity is of size(" << intVal_.size() << ")");
  intVal_.erase(intVal_.begin() + ipart);
}

void Atom::swap(const int previous, const int finl) {
  const int tmp = intVal_[finl];
  intVal_[finl] = intVal_[previous];
  intVal_[previous] = tmp;
}

void Atom::initDataType(const std::string type) {
  ASSERT(type == "int", "only implemented for integer types.");
  type_ = type;
}

unsigned int Atom::size() const {
  return intVal_.size();
}

void Atom::writeRestart(const char* fileName) {
  fileBackUp(fileName);
  std::ofstream file(fileName);
  file << "# datatype " << type_ << endl;
  vecRestartPrinter("intVal", intVal_, fileName);
}

Atom::Atom(const char* fileName) {
  ASSERT(fileExists(fileName), "restart file(" << fileName
         << ") doesn't exist");
  defaultConstruction_();
  type_ = fstos("datatype", fileName);
  vecRestartReader("intVal", &intVal_, fileName);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

