
#ifndef FEASST_CHAIN_SELECT_REPTATE_H_
#define FEASST_CHAIN_SELECT_REPTATE_H_

#include "configuration/include/select.h"
#include "chain/include/select_end_segment.h"

namespace feasst {

/// Select a random end point for reptation.
class SelectReptate : public SelectEndSegment {
 public:
  /**
    args:
    - SelectEndSegment arguments.
   */
  explicit SelectReptate(argtype args = argtype());
  explicit SelectReptate(argtype * args);

  void precompute(System * system) override;

  void update_anchor(const bool is_endpoint_beginning,
    const System * system) override;

  void mid_stage() override;

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
    argtype args = argtype()) {
  return std::make_shared<SelectReptate>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_SELECT_REPTATE_H_
