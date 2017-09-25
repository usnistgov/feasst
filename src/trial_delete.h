#ifndef TRIAL_DELETE_H_
#define TRIAL_DELETE_H_

#include <memory>
#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class Space;
class Pair;
class Criteria;

/**
 * Attempt to delete particle(s).
 */
class TrialDelete : public Trial {
 public:
  /// Constructor
  TrialDelete(Space *space, Pair *pair, Criteria *criteria, const char* molType);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  explicit TrialDelete(const char* molType);

  /// Constructors with no molType listed will attempt to delete any particle.
  TrialDelete(Space *space, Pair *pair, Criteria *criteria);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  TrialDelete();

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  TrialDelete(const char* fileName, Space *space, Pair *pair,
              Criteria *criteria);

  ~TrialDelete() {}
  TrialDelete* clone(Space* space, Pair *pair, Criteria *criteria) const {
    TrialDelete* t = new TrialDelete(*this);
    t->reconstruct(space, pair, criteria); return t; }
  shared_ptr<TrialDelete> cloneShrPtr
    (Space* space, Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialDelete, Trial>
    (cloneImpl(space, pair, criteria))); }

 protected:
  string molType_;   //!< type of molecule to delete
  int molid_;        //!< index of molecule type

  /// attempt deletion
  void attempt1_();

  // default constructor
  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Space* space, Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialDelete> t = make_shared<TrialDelete>(*this);
    t->reconstruct(space, pair, criteria); return t;
  }
};

class MC;

/// Add a "TrialDelete" object to the Monte Carlo object, mc
void deleteTrial(MC *mc, const char* moltype);

/// Add a "TrialDelete" object to the Monte Carlo object, mc
void deleteTrial(shared_ptr<MC> mc, const char* moltype);

/// Add a "TrialDelete" object to the Monte Carlo object, mc
void deleteTrial(MC *mc);

/// Add a "TrialDelete" object to the Monte Carlo object, mc
void deleteTrial(shared_ptr<MC> mc);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_DELETE_H_

