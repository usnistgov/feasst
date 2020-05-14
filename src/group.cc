/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./group.h"
#include "./space.h"

namespace feasst {

Group::Group() : Atom() {
  initDataType("int");
}

void Group::molid(const int iMol) {
  molid_.push_back(iMol);
}

//void Group::moltype(const char* name) {
//  moltype_.push_back(name);
//}

int Group::iatomInGroup_(const Space &space, const int iatom) {
  const int iMol = space.mol()[iatom];
  for (vector<int>::iterator it = molid_.begin();
       it != molid_.end();
       ++it) {
    if (space.molid()[iMol] == *it) {
      return 1;
    }
  }
  const std::string molType = space.moltype()[iMol];
  for (vector<std::string>::iterator it = moltype_.begin();
       it != moltype_.end();
       ++it) {
    if (molType == *it) {
      return 1;
    }
  }
  return 0;
}

void Group::writeRestart(const char* fileName) {
  Atom::writeRestart(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  vecRestartPrinter("molid", molid_, fileName);
}

Group::Group(const char* fileName)
  : Atom(fileName) {
  vecRestartReader("molid", &molid_, fileName);
  // vecRestartReader("moltype", &moltype_, fileName);
}

void Group::addPart(const Space &space) {
//  const int nAtomPrev = intVal_.size();
  Atom::addPart(space);

//  // find the last mol that was added via addMol(""), add to subSpace
//  const int molid = space->molid()[nAtomPrev];
//  const std::string moltype = space->moltype()[molid];
//  subSpace_->addMol(moltype);
//
//  // update subspace and per-atom-quantity
//  if (type_ == "int") {
//    for (int iatom = space->natom() - 1; iatom >= nAtomPrev; --iatom) {
//      if (intVal_[iatom] == 0) {
//        subSpace_->delPart(iatom);
//      } else {
//        atom.setVal(iatom, iatom);
//      }
//    }
//  }
}

void Group::delPart(const int ipart) {
  Atom::delPart(ipart);
  // update based on new per-atom-quantity
  // broken because can't know how molecules will be treated on atom basis
}

void Group::initSubSpace(shared_ptr<Space> space) {
  subSpace_ = space;

  // add a new per-atom-quantity to space, which relates the id's of the
  // particles in original space to the subspace
  auto atom = make_shared<Atom>();
  atom->initName("subspace");
  atom->initDataType("int");
  subSpace_->initAtom(atom);

  // initialize particles that are in space
  addPart(*subSpace_);

  // remove atoms in subspace which are not part of the group
  // otherwise, label the per-atom-quantity for the full space index
  const int nAtomPrev = intVal_.size();
  if (type_ == "int") {
    for (int iatom = subSpace_->natom() - 1; iatom >= nAtomPrev; --iatom) {
      if (intVal_[iatom] == 0) {
        subSpace_->delPart(iatom);
      } else {
        atom->setVal(iatom, iatom);
      }
    }
  }
}

vector<int> Group::randMol(Space * space) {
  // to begin, select a few random particles. If those are not in the group,
  // then attempt a different method
  vector<int> mpart;
  int nAttempts = 0;
  int found = 0;
  while ((found == 0) && (nAttempts < 10)) {
    mpart = space->randMol();
    if (intVal_[mpart[0]] == 1) {
      found = 1;
    }
  }

  if (found == 0) {
    shared_ptr<Space> spaceShrPtr = make_shared<Space>(*space);
    initSubSpace(spaceShrPtr);
    mpart = subSpace_->randMol();
  }
  return mpart;
}

}  // namespace feasst

