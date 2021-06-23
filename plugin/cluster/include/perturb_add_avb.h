
#ifndef FEASST_CLUSTER_PERTURB_ADD_AVB_H_
#define FEASST_CLUSTER_PERTURB_ADD_AVB_H_

#include "cluster/include/perturb_move_avb.h"

namespace feasst {

// HWH optimize -> update cell list in finalize only?
/**
  Add a particle within the AV of an existing particle.
 */
class PerturbAddAVB : public Perturb {
 public:
  /**
    args:
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
    - delay_add: If true, don't add particle until finalize (default: true).
   */
  PerturbAddAVB(argtype args = argtype());
  PerturbAddAVB(argtype * args);

  /// Return if the particle isn't added until finalized.
  bool delay_add() const { return delay_add_; }

  void precompute(TrialSelect * select, System * system) override {
    select->set_ghost(true); }

  void perturb(
      System * system,
      TrialSelect * select,
      Random * random,
      const bool is_position_held = false) override;

  void revert(System * system) override;
  void finalize(System * system) override;
  std::string status_header() const override;
  std::string status() const override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbAddAVB(std::istream& istr);
  virtual ~PerturbAddAVB() {}

 private:
  bool delay_add_;
  std::shared_ptr<PerturbMoveAVB> move_;
};

inline std::shared_ptr<PerturbAddAVB> MakePerturbAddAVB(
    argtype args = argtype()) {
  return std::make_shared<PerturbAddAVB>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_PERTURB_ADD_AVB_H_
