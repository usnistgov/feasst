
#ifndef FEASST_CHAIN_SELECT_REPTATE_H_
#define FEASST_CHAIN_SELECT_REPTATE_H_

#include "chain/include/select_end_segment.h"

namespace feasst {

/// Select a random end point for reptation.
class SelectReptate : public SelectEndSegment {
 public:
  SelectReptate(const argtype& args = argtype());

  void precompute(System * system) override;

  void update_anchor(const bool is_endpoint_beginning,
    const System * system) override;

  void mid_stage() override {
    // exclude the anchor from interactions.
    // include interactions with site that use to be bonded
    mobile_.set_new_bond(anchor_);
    mobile_.set_old_bond(bonded_to_);
  }

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectReptate(std::istream& istr);
  virtual ~SelectReptate() {}

 protected:
  void serialize_select_reptate_(std::ostream& ostr) const;

 private:
  Select bonded_to_;
};

inline std::shared_ptr<SelectReptate> MakeSelectReptate(
    const argtype& args = argtype()) {
  return std::make_shared<SelectReptate>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_SELECT_REPTATE_H_
