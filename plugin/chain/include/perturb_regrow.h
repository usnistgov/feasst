
#ifndef FEASST_CHAIN_PERTURB_REGROW_H_
#define FEASST_CHAIN_PERTURB_REGROW_H_

#include "core/include/perturb_move.h"

namespace feasst {

/**
 */
class PerturbRegrow : public PerturbSelectMove {
 public:
  // assumes that if selection begins with first site, then start regrow
  // from opposite end (by using custom flag)
  // assumes linear chain
  void perturb(System * system) override {
    Configuration * config = get_config_before_move(system);
    SelectList bonded = selection();
    DEBUG("regrowing sites " << selection().str());
    ASSERT(selection().num_sites() == anchor().num_sites(), "error");
//    DEBUG("is_reversed " << is_custom_flag());
//    const int num_sites = static_cast<int>(bonded.site_indices()[0].size());
    ASSERT(selection().num_sites() == 1, "error");
//    for (int index = 0; index < selection().num_sites(); ++index) {
//      int to_update = index;
//      if (is_custom_flag()) {
//        to_update = selection().num_sites() - 1 - index;
//      }
      rebond_(0, system, &bonded);
//    }

    // recenter particle position
    bonded.set_particle_position(0,
      system->configuration().select_particle(0).average_site_position());

    config->update_positions(bonded);
    after_move();
  }

  void set_in_sphere(const double bond_length, const int index, SelectList * bonded);

  ~PerturbRegrow() {}

 private:
  Random random_;

  void rebond_(const int site_to_update,
               System * system,
               SelectList * bonded);
};

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_REGROW_H_
