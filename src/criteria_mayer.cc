#include "./criteria_mayer.h"
#include "./functions.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

CriteriaMayer::CriteriaMayer(const double beta
  ) : Criteria(beta, 0.) {
  defaultConstruction_();
}

CriteriaMayer::CriteriaMayer(const char* fileName)
  : Criteria(fileName) {
  defaultConstruction_();
}

void CriteriaMayer::defaultConstruction_() {
  className_.assign("CriteriaMayer");
  verbose_ = 0;
}

CriteriaMayer* CriteriaMayer::clone() const {
  CriteriaMayer* c = new CriteriaMayer(*this);
  c->reconstruct();
  return c;
}

shared_ptr<CriteriaMayer> CriteriaMayer::cloneShrPtr() const {
  return(std::static_pointer_cast<CriteriaMayer, Criteria>(cloneImpl_()));
}

shared_ptr<Criteria> CriteriaMayer::cloneImpl_() const {
  shared_ptr<CriteriaMayer> c = make_shared<CriteriaMayer>(*this);
  c->reconstruct();
  return c;
}

int CriteriaMayer::accept(const double lnpMet, const double peNew,
  const char* moveType, const int reject) {
  string moveTypeStr(moveType);
  ASSERT(moveTypeStr == "move", "unrecognized moveType("
    << moveType << ")");
  
  // acceptance criteria
  const double f12 = exp(-beta_*peNew) - 1.,
    f12old = exp(-beta_*peOld_) - 1.;
//  cout << "lnpMEt " << lnpMet << " peold " << peOld_ << " peNew " << peNew << " penewcalc " << exp(lnpMet)/-beta_ << endl;
  int returnVal = -1;
  double f12current, f12ref;
  if ( (reject != 1) && 
       (uniformRanNum() < fabs(f12)/fabs(f12old)) ) {
    returnVal = 1;
    f12current = f12;
    pairRef_->initEnergy();
    f12ref = exp(-beta_*pairRef_->peTot()) - 1.;
  } else {
    returnVal = 0;
    f12current = f12old;
    f12ref = exp(-beta_*peRefOld_) - 1.;
  }

  // accumulate mayer functions
  if (f12current < 0) {
    mayer.accumulate(-1);
  } else {
    mayer.accumulate(1);
  }
  mayerRef.accumulate(f12ref/fabs(f12current));

//  cout << f12current << " " << f12ref << " " << mayer.average() << " "
//    << mayerRef.average() << " " << mayer.average()/mayerRef.average() << " "
//    << "returnVal " << returnVal << " prob " << fabs(f12)/fabs(f12old) << " "
//    << "peNew " << peNew << endl;

  // return acceptance
  ASSERT( (returnVal == 0) || (returnVal == 1), "returnVal(" << returnVal
    << "not 0 or 1");
  return returnVal;
}

void CriteriaMayer::store(const Space* space, Pair* pair) {
  if (space == NULL) {}
  peOld_ = pair->peTot();
  ASSERT(fabs(peOld_) > 0, "Mayer sampling requires that the particles never"
    << "leave their potential wells, e.g., the potential energy is nonzero");
  pairRef_->initEnergy();
  peRefOld_ = pairRef_->peTot();
  ASSERT(space == pairRef_->space(), "reference potential must point to the"
    << "same space object");
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_


