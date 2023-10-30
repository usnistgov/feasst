
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
#include "aniso/include/rotator.h"

namespace feasst {

/**
  Generate a table of interactions between two rigid bodies in 3D.
  This class is currently in development and not well supported.

  Contact distances are determined by RefPotential.
  This is typically a HardSphere potential with an optimzied VisitModelCell.
 */
class TabulateTwoRigidBody3D : public Action {
 public:
  //@{
  /** @name Arguments
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

  //@}
  /** @name Public Functions
   */
  //@{

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
  //@}
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
