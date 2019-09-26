
#ifndef FEASST_CHAIN_PERTURB_REPTATE_H_
#define FEASST_CHAIN_PERTURB_REPTATE_H_

#include "monte_carlo/include/perturb_distance.h"

namespace feasst {

/**
  For a reptation, if new bond is accepted, then change the positions of all the
  sites along the chain.
 */
class PerturbReptate : public PerturbDistance {
 public:
  PerturbReptate(const argtype& args = argtype()) : PerturbDistance(args) {
    class_name_ = "PerturbReptate";
  }
  void move(System * system, TrialSelect * select, Random * random) override {
    PerturbDistance::move(system, select, random);
    set_finalize_possible(true, select);
  }

  void finalize(System * system) override {
    const SelectList& mobile = finalize_select()->mobile();
    const int part_index = mobile.particle_indices()[0];
    SelectList entire = SelectList().particle(part_index,
                                              system->configuration(),
                                              0 // group that includes all
                                              );
    if (mobile.site_indices()[0][0] == 0) {
      for (int site = 1; site < entire.num_sites(); ++site) {
        entire.set_site_position(0, site - 1, entire.site_positions()[0][site]);
        entire.set_site_properties(0, site - 1, entire.site_properties()[0][site]);
      }
      entire.set_site_position(0, entire.num_sites() - 1, mobile.site_positions()[0][0]);
      entire.set_site_properties(0, entire.num_sites() - 1, mobile.site_properties()[0][0]);
    } else {
      for (int site = entire.num_sites() - 1; site >= 1; --site) {
        entire.set_site_position(0, site, entire.site_positions()[0][site - 1]);
        entire.set_site_properties(0, site, entire.site_properties()[0][site - 1]);
      }
      entire.set_site_position(0, 0, mobile.site_positions()[0][0]);
      entire.set_site_properties(0, 0, mobile.site_properties()[0][0]);
    }
    DEBUG("entire " << entire.str() << " pos " << entire.site_positions()[0][0].str() << " end " << entire.site_positions()[0][49].str());
    system->get_configuration()->update_positions(entire, false);
  }
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbReptate(std::istream& istr);
  virtual ~PerturbReptate() {}

 protected:
  void serialize_perturb_reptate_(std::ostream& ostr) const;
};

inline std::shared_ptr<PerturbReptate> MakePerturbReptate(const argtype& args = argtype()) {
  return std::make_shared<PerturbReptate>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_REPTATE_H_
