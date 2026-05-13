
#ifndef FEASST_ANISO_ROTATOR_H_
#define FEASST_ANISO_ROTATOR_H_

#include <vector>
#include <memory>
#include "math/include/formula.h"
#include "configuration/include/file_xyz.h"
#include "system/include/system.h"
#include "monte_carlo/include/action.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "monte_carlo/include/perturb_translate.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class Rotator {
 public:
  //@{
  /** @name Arguments
    - unique_tolerance: tolerance to determine uniqueness based on positions
      (default: 1e-5);
    - contact_tolerance: contact distance tolerance (default: 1e-4).
    - num_proc: number of processors (default: 1).
    - proc: index of processor (default: 0).
    - hard_limit_u: potential energy that is considered contact.
      A user may set to positive values.
      By default, it is set to 1e10, which uses an objective function
      specificly for the HardSphere Potential.
   */
  explicit Rotator(argtype args = argtype());
  explicit Rotator(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

//  bool identical(const Configuration& config) const;
  std::vector<std::vector<double> > gen_global_bounds(const Configuration& config) const;
  void gen_orientations(const int num_orientations_per_pi, const Configuration& config,
    /// optional custom bounds. If empty, use full bounds.
    std::vector<std::vector<double> > bounds = std::vector<std::vector<double> >());
  void gen_unique_orientations(const int num_orientations_per_pi, System * system,
    std::vector<std::vector<double> > bounds = std::vector<std::vector<double> >());
  void init(System * system, const std::string xyz_file,
    const std::string contact_xyz_file);
  int num_orientations() const { return static_cast<int>(stheta_.size()); }
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
  double energy(const int ior, const double displacement, System * system);
  double contact_distance(const int ior, System * system);
  double cutoff_distance(const int ior, System * system);
  double hard_limit_u_;
  constexpr static double hard_u_ = 1e10; // ensure this matches hard_limit_u doc
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
  std::vector<float> cutoff_;
  std::vector<std::vector<float> > energy_;
  std::vector<std::vector<int> > indices_;
  std::vector<int> sizes_;
  int num_proc_;
  int proc_;
  int num_orientations_all_proc_;

  // uniqueness
  std::vector<std::vector<Position> > last_three_sites_;
  std::vector<Position> last_three_;
  std::vector<double> tmp_vec_;
  std::vector<int> unique_;
  double unique_tolerance_;
  double contact_tolerance_;
  Position tmp1_, tmp2_;
  //@}
};

}  // namespace feasst

#endif  // FEASST_ANISO_ROTATOR_H_
