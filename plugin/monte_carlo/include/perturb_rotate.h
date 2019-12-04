
#ifndef FEASST_MONTE_CARLO_PERTURB_ROTATE_H_
#define FEASST_MONTE_CARLO_PERTURB_ROTATE_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

/**
  Rotate the positions of the selection.
 */
class PerturbRotate : public PerturbMove {
 public:
  PerturbRotate(const argtype& args = argtype());

  /// Change the position in the selection given a pivot and rotation matrix.
  void update_selection(const Position& pivot,
      const RotationMatrix& rotation,
      TrialSelect * select,
      /// Rotate particle positions (default). Otherwise, do not.
      const bool rotate_particle_position = true);

  /// Change the position of the selection given a pivot and rotation matrix.
  void move(const Position& pivot,
      const RotationMatrix& rotation,
      System * system,
      TrialSelect * select,
      /// Rotate particle positions (default). Otherwise, do not.
      const bool rotate_particle_position = true);

  /// Rotate the selected particles using the tuning parameter.
  /// Set the pivot using the first particle position, and also
  /// rotate the particle positions.
  void move(System * system,
      TrialSelect * select,
      Random * random) override;

  /// Rotate the selected particles using the tuning parameter.
  void move(System * system,
      TrialSelect * select,
      Random * random,
      /// If pivot is empty, use first particle position.
      const Position& pivot,
      /// Rotate particle positions if true. Otherwise, do not.
      const bool rotate_particle_position);

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbRotate(std::istream& istr);
  virtual ~PerturbRotate() {}

 protected:
  void serialize_perturb_rotate_(std::ostream& ostr) const;

 private:
  const Position& piv_sel_(const Position& pivot, const TrialSelect * select) {
    if (pivot.dimension() == 0) {
      return select->mobile().particle_positions()[0];
    }
    return pivot;
  }

  // optimization
  // check to see if rotation is not necessary when pivot is equivalent to
  // the only rotated position.
  bool is_rotation_not_needed_(const TrialSelect * select,
      const Position& pivot) {
    const SelectList& rotated = select->mobile();
    if (rotated.num_particles() == 1) {
      if (static_cast<int>(rotated.site_indices()[0].size()) == 1) {
        if (rotated.site_positions()[0][0].is_equal(pivot)) {
          return true;
        }
      }
    }
    return false;
  }
};

inline std::shared_ptr<PerturbRotate> MakePerturbRotate(
    const argtype& args = argtype()) {
  return std::make_shared<PerturbRotate>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_ROTATE_H_
