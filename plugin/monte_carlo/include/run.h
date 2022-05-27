
#ifndef FEASST_MONTE_CARLO_RUN_H_
#define FEASST_MONTE_CARLO_RUN_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Perform a number of trials.
 */
class Run : public Action {
 public:
  /**
    The following arguments are completed in the order listed.

    args:
    - num_trials: run this many trials (default: -1. e.g., None)
    - until_num_particles: run until this many particles (default: -1. e.g., None)
    - particle_type: type of particle to count. If -1, all particles (default: -1).
    - for_hours: run for this many CPU hours (default: -1 e.g., None).
    - until_criteria_complete: run until Criteria is complete (default: false)
   */
  explicit Run(argtype args = argtype());
  explicit Run(argtype * args);
  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<Run>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<Run>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Run(std::istream& istr);
  virtual ~Run() {}

 private:
  int num_trials_;
  int until_num_particles_;
  int particle_type_;
  double for_hours_;
  bool until_criteria_complete_;
};

inline std::shared_ptr<Run> MakeRun(argtype args = argtype()) {
  return std::make_shared<Run>(args);
}

/**
  Remove a Trial.
 */
class RemoveTrial : public Action {
 public:
  /**
    args:
    - index: index of trial to remove, in order of trials added.
      If -1, do nothing. (default: -1).
    - name: remove first trial with this class name, if not empty.
      (default: empty).
    - all: remove all trials (default: false)
   */
  explicit RemoveTrial(argtype args = argtype());
  explicit RemoveTrial(argtype * args);
  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<RemoveTrial>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<RemoveTrial>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RemoveTrial(std::istream& istr);
  virtual ~RemoveTrial() {}

 private:
  int index_;
  std::string name_;
  bool all_;
};

inline std::shared_ptr<RemoveTrial> MakeRemoveTrial(argtype args = argtype()) {
  return std::make_shared<RemoveTrial>(args);
}

/**
  Remove a Analyze.
 */
class RemoveAnalyze : public Action {
 public:
  /**
    args:
    - index: index of analyze to remove, in the order added.
      If -1, do nothing. (default: -1).
    - name: remove first analyze with this class name, if not empty.
      (default: empty).
    - all: remove all analyzers (default: false)
   */
  explicit RemoveAnalyze(argtype args = argtype());
  explicit RemoveAnalyze(argtype * args);
  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<RemoveAnalyze>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<RemoveAnalyze>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RemoveAnalyze(std::istream& istr);
  virtual ~RemoveAnalyze() {}

 private:
  int index_;
  std::string name_;
  bool all_;
};

inline std::shared_ptr<RemoveAnalyze> MakeRemoveAnalyze(argtype args = argtype()) {
  return std::make_shared<RemoveAnalyze>(args);
}

/**
  Remove a Modify.
 */
class RemoveModify : public Action {
 public:
  /**
    args:
    - index: index of modify to remove, in the order added.
      If -1, do nothing. (default: -1).
    - name: remove first modify with this class name, if not empty.
      (default: empty).
    - all: remove all modifiers (default: false)
   */
  explicit RemoveModify(argtype args = argtype());
  explicit RemoveModify(argtype * args);
  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<RemoveModify>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<RemoveModify>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RemoveModify(std::istream& istr);
  virtual ~RemoveModify() {}

 private:
  int index_;
  std::string name_;
  bool all_;
};

inline std::shared_ptr<RemoveModify> MakeRemoveModify(argtype args = argtype()) {
  return std::make_shared<RemoveModify>(args);
}

/**
  Write a Checkpoint.
 */
class WriteCheckpoint : public Action {
 public:
  /**
    args:
   */
  explicit WriteCheckpoint(argtype args = argtype());
  explicit WriteCheckpoint(argtype * args);
  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<WriteCheckpoint>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<WriteCheckpoint>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit WriteCheckpoint(std::istream& istr);
  virtual ~WriteCheckpoint() {}
};

inline std::shared_ptr<WriteCheckpoint> MakeWriteCheckpoint(argtype args = argtype()) {
  return std::make_shared<WriteCheckpoint>(args);
}

/**
  Make a new reference potential based on an existing potential.
 */
class ConvertToRefPotential : public Action {
 public:
  /**
    args:
    - potential_index: index of the full potential to copy as a template
      (default: 0).
    - cutoff: set cutoff of all site_types to this value.
      Ignore if -1 (default: -1).
    - use_cell: use VisitModelCell with cutoff as min_length (default: false).
   */
  explicit ConvertToRefPotential(argtype args = argtype());
  explicit ConvertToRefPotential(argtype * args);
  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<ConvertToRefPotential>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<ConvertToRefPotential>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ConvertToRefPotential(std::istream& istr);
  virtual ~ConvertToRefPotential() {}

 private:
  int potential_index_;
  double cutoff_;
  bool use_cell_;
};

inline std::shared_ptr<ConvertToRefPotential> MakeConvertToRefPotential(
    argtype args = argtype()) {
  return std::make_shared<ConvertToRefPotential>(args);
}

/**
  Add a reference potential.
 */
class RefPotential : public Action {
 public:
  /**
    args:
    - Same arguments as Potential.
   */
  explicit RefPotential(argtype args = argtype());
  explicit RefPotential(argtype * args);
  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<RefPotential>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<RefPotential>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RefPotential(std::istream& istr);
  virtual ~RefPotential() {}

 private:
  argtype args_;
};

inline std::shared_ptr<RefPotential> MakeRefPotential(
    argtype args = argtype()) {
  return std::make_shared<RefPotential>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_RUN_H_
