
#ifndef FEASST_MONTE_CARLO_PERTURB_ADD_H_
#define FEASST_MONTE_CARLO_PERTURB_ADD_H_

#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

// HWH optimize -> update cell list in finalize only?
/**
  Add a particle to the system.
 */
class PerturbAdd : public Perturb {
 public:
  /**
    args:
    - delay_add: If true, don't add particle until finalize (default: true).
   */
  explicit PerturbAdd(argtype args = argtype());
  explicit PerturbAdd(argtype * args);

  /// Return if the particle isn't added until finalized.
  bool delay_add() const { return delay_add_; }

  //initialize ghost selection in TrialSelect?
  void precompute(TrialSelect * select, System * system) override;

  void before_select() override;

  void perturb(
      System * system,
      TrialSelect * select,
      Random * random,
      const bool is_position_held = false,
      Acceptance * acceptance = NULL) override {
    add(system, select, random, empty_, is_position_held); }

  /// Add select to the system.
  void add(
    System * system,
    TrialSelect * select,
    Random * random,
    /// place particle anywhere if center is of zero dimension.
    const Position& center,
    const bool is_position_held = false);

  void revert(System * system) override;
  void finalize(System * system) override;
  std::string status_header() const override;
  std::string status() const override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbAdd(std::istream& istr);
  virtual ~PerturbAdd() {}

 private:
  bool delay_add_;

  // temporary
  Position empty_;
  PerturbAnywhere anywhere_;
};

inline std::shared_ptr<PerturbAdd> MakePerturbAdd(argtype args = argtype()) {
  return std::make_shared<PerturbAdd>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_ADD_H_
