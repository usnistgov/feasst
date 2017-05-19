#include "./pair_wall.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairWall::PairWall(Space* space, Barrier* barrier)
  : Pair(space, 0.),
    barrier_(barrier) {
  className_.assign("PairWall");
}

double PairWall::multiPartEner(const vector<int> mpart, const int flag) {
  if (flag == 0) {}; //remove unused parameter warning
  const vector<double> &x = space_->x();
  peSRone_ = 0;
  for (int impart = 0; impart < int(mpart.size()); ++impart) {
    const int iPart = mpart[impart];
    const int iType = space_->type()[iPart];

    vector<double> coord(x.begin() + iPart*dimen_, 
                         x.begin() + (iPart+1)*dimen_);
    peSRone_ += barrier_->potential(coord, sig_[iType]);
    
    // here is where the barrier class would be implemented
    // one possibility is to take what is in Mahynski's "system" class
    // and put in "space", then access the space variables space_->
    // or implement a container for multiple barriers
  }
  return peSRone_;
}

int PairWall::initEnergy() {
  std::fill(pe_.begin(), pe_.end(), 0.);
  fill(0., f_);
  fill(0., vr_);
  // compute energy
  multiPartEner(space_->listAtoms());
  peTot_  = peSRone_;
  return 0;
}

void PairWall::update(const vector<int> mpart,    //!< particles involved in move
                    const int flag,         //!< type of move
                    const char* uptype    //!< description of update type
  ) {
  std::string uptypestr(uptype);
  if (uptypestr.compare("store") == 0) {
    if (flag == 0 || flag == 2 || flag == 3) {
      deWall_ = peSRone_;
    }
  }

  if (uptypestr.compare("update") == 0) {
    if (flag == 0) {
      peWall_ += peSRone_ - deWall_;
      peTot_ += peSRone_ - deWall_;
    }
    if (flag == 2) {
      peWall_ -= deWall_;
      peTot_ -= deWall_;
    }
    if (flag == 3) {
      peWall_ += deWall_;
      peTot_ += deWall_;
    }
  }
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_
