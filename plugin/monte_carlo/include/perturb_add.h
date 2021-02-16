
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
  explicit PerturbAdd(argtype args = argtype());
  explicit PerturbAdd(argtype * args);

  //initialize ghost selection in TrialSelect?
  void precompute(TrialSelect * select, System * system) override {
    select->set_ghost(true); }

  void perturb(
      System * system,
      TrialSelect * select,
      Random * random,
      const bool is_position_held = false) override {
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
  // temporary
  Position empty_;
  PerturbAnywhere anywhere_;
};

inline std::shared_ptr<PerturbAdd> MakePerturbAdd(argtype args = argtype()) {
  return std::make_shared<PerturbAdd>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_ADD_H_
