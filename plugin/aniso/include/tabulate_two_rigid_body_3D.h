
#ifndef FEASST_ANISO_TABULATE_TWO_RIGID_BODY_3D_H_
#define FEASST_ANISO_TABULATE_TWO_RIGID_BODY_3D_H_

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
  /**
    args:
    - unique_tolerance: tolerance to determine uniqueness based on positions
      (default: 1e-5);
    - contact_tolerance: contact distance tolerance (default: 1e-4).
    - num_proc: number of processors (default: 1).
    - proc: index of processor (default: 0).
   */
  explicit Rotator(argtype args = argtype());
  explicit Rotator(argtype * args);
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
};

/**
  Generate a table of interactions between two rigid bodies in 3D.
  This class is currently in development and not well supported.

  Contact distances are determined by RefPotential.
  This is typically a HardSphere potential with an optimzied VisitModelCell.
 */
class TabulateTwoRigidBody3D : public Action {
 public:
  /**
    args:
    - num_orientations_per_pi: As described in VisitModelInnerTable
      (default: -1).
    - num_z: As described in VisitModelInnerTable (default: -1).
      Set to -1 when computing contact distances.
    - gamma: As described in VisitModelInnerTable (default: -4).
    - smoothing_distance: As described in VisitModelInnerTable (default: 2).
    - max_energy: For handshake configurations (large energy for z > 0),
      set the energy to max_energy_set instead, which may be 1e30, near the max
      single precision, or lower (default: 1e30).
    - max_energy_set: set to this energy if above max and z > 0 (default: 5).
    - output_orientation_file: name of file to output unique orientations,
      and terminate early, if not empty (default: empty).
      This is used to compute unique orientations using a fast and simple model.
    - input_orientation_file: name of file to input unique orientations,
      if not empty (default: empty).
    - output_table_file: name of file to output computed table.
    - input_table_file: name of file to input computed table.
    - xyz_file: if not empty, visualize (default: empty).
    - contact_xyz_file: if not empty, visualize (default: empty).
    - Rotator arguments.
   */
  explicit TabulateTwoRigidBody3D(argtype args = argtype());
  explicit TabulateTwoRigidBody3D(argtype * args);

  /// Return the maximum cubic side length required for particle type.
  double max_cubic_side_length(const int particle_type,
                               const Configuration& config) const;

  /// Adjust the Domain.
  void adjust_domain(System * system);

  const Rotator& rotator() const { return rotator_; }

  /**
    Write the table file as described in VisitModelInnerTable.
    If num_z == -1, include the unique id for use in subsequent energy calc.
   */
  void write_table(MonteCarlo * mc) const;

  /**
    Read table for contact (num_z == -1) produced from the above,
    in order to compute energy table (num_z > 1).
   */
  void read_contact_table(MonteCarlo * mc);

  /**
    Write a combined table file.

    args:
    - prefix: characters in file name before processor index.
    - suffix: characters in file name after processor index.
   */
  void combine_table(argtype args = argtype()) const;

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<TabulateTwoRigidBody3D>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<TabulateTwoRigidBody3D>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TabulateTwoRigidBody3D(std::istream& istr);
  virtual ~TabulateTwoRigidBody3D() {}

  void read_input_orientations_(const Configuration& config);
 private:
  int num_orientations_per_pi_;
  int num_z_;
  double gamma_;
  double smoothing_distance_;
  double max_energy_;
  double max_energy_set_;
  std::string output_orientation_file_;
  std::string input_orientation_file_;
  std::string output_table_file_;
  std::string input_table_file_;
  std::string xyz_file_;
  std::string contact_xyz_file_;
  Rotator rotator_;
  void ouput_orientations_();
};

inline std::shared_ptr<TabulateTwoRigidBody3D> MakeTabulateTwoRigidBody3D(argtype args = argtype()) {
  return std::make_shared<TabulateTwoRigidBody3D>(args);
}

}  // namespace feasst

#endif  // FEASST_ANISO_TABULATE_TWO_RIGID_BODY_3D_H_
