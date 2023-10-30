
#ifndef FEASST_ANISO_ROTATOR_H_
#define FEASST_ANISO_ROTATOR_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "math/include/formula.h"
#include "configuration/include/file_xyz.h"
#include "system/include/system.h"
#include "monte_carlo/include/action.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "monte_carlo/include/perturb_translate.h"

namespace feasst {

class Rotator;

class ContactObjective : public Formula {
 public:
  ContactObjective(Rotator * rotator, System * system, const double ior, const int ref_pot);
  double evaluate(const double distance) const override;
 private:
  Rotator * rotator_;
  System * system_;
  int ior_;
  int ref_pot_;
};

class Rotator {
 public:
  //@{
  /** @name Arguments
    - unique_tolerance: tolerance to determine uniqueness based on positions
      (default: 1e-5);
    - contact_tolerance: contact distance tolerance (default: 1e-4).
    - num_proc: number of processors (default: 1).
    - proc: index of processor (default: 0).
   */
  explicit Rotator(argtype args = argtype());
  explicit Rotator(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void gen_orientations(const int num_orientations_per_pi, const Configuration& config);
  void init(System * system, const std::string xyz_file,
    const std::string contact_xyz_file);
  int num_orientations() const { return static_cast<int>(eulers_.size()); }
  int num_orientations_all_proc() const { return num_orientations_all_proc_; }
  void determine_if_unique(const int ior, const std::vector<int>& iors, const int num_threads, System * system);
  void update_xyz(const int ior, const double displacement, System * system);
  void set_last_three_sites(const int ior, System * system);
  void check_last_three_sites(const int ior, System * system);
  void revert(System * system);
  int num_unique() const;
  int num_proc() const { return num_proc_; }
  bool ior_in_proc(const int ior) const;
//  int ior_all_proc_to_ior(const int iorall) const {
//    if (ior_in_proc(iorall)) {
//      return int(iorall/num_proc_);
//    }
//    return -1;
//  }
  double fraction_unique() const {
    return num_unique()/static_cast<double>(unique_.size()); };
  double energy(const int ior, const double displacement, System * system, const int ref_potential = -1);
  double contact_distance(const int ior, System * system, const int ref_potential = 0);
  std::string xyz_file_name_;
  std::string contact_xyz_file_name_;
  Position com1_;
  std::shared_ptr<TrialSelectParticle> select_;
  std::shared_ptr<PerturbRotate> rotate_;
  std::shared_ptr<PerturbTranslate> translate_;
  std::shared_ptr<Position> origin_;
  RotationMatrix rot_mat_;
  FileXYZ xyz_;
  FileXYZ contact_f_;
  std::vector<double> stheta_;
  std::vector<double> sphi_;
  std::vector<Euler> eulers_;
  std::vector<float> contact_;
  std::vector<std::vector<float> > energy_;
  int num_proc_;
  int proc_;
  int num_orientations_all_proc_;

  // uniqueness
  std::vector<std::vector<Position> > last_three_sites_;
  std::vector<Position> last_three_;
  std::vector<double> tmp_3vec_;
  std::vector<int> unique_;
  double unique_tolerance_;
  double contact_tolerance_;
  Position tmp1_, tmp2_;
  //@}
};

}  // namespace feasst

#endif  // FEASST_ANISO_ROTATOR_H_
