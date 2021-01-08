
#ifndef FEASST_MONTE_CARLO_PERTURB_ROTATE_H_
#define FEASST_MONTE_CARLO_PERTURB_ROTATE_H_

#include "math/include/matrix.h"
#include "monte_carlo/include/perturb_move.h"

namespace feasst {

/**
  Rotate the positions of the selection.
  Assumes the selection is only one particle.
 */
class PerturbRotate : public PerturbMove {
 public:
  /**
    args:
    - pivot_site: set the site index in selection with which to use as the
      pivot for rotation (default: 0).
   */
  PerturbRotate(const argtype& args = argtype());

  /// Change the position in the selection given a pivot and rotation matrix.
  void update_selection(const Position& pivot,
      const RotationMatrix& rotation,
      TrialSelect * select);

  /// Change the position of the selection given a pivot and rotation matrix.
  void move(const Position& pivot,
      const RotationMatrix& rotation,
      System * system,
      TrialSelect * select);

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
      const Position& pivot);

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbRotate(std::istream& istr);
  virtual ~PerturbRotate() {}

 protected:
  void serialize_perturb_rotate_(std::ostream& ostr) const;

 private:
  int pivot_site_;

  // temporary and not serialized
  RotationMatrix rot_mat_tmp_;
  Position axis_tmp_;

  const Position& piv_sel_(const Position& pivot, const TrialSelect * select) {
    if (pivot.dimension() == 0) {
      return select->mobile().site_positions()[0][0];
    }
    return pivot;
  }

  // optimization
  // check to see if rotation is not necessary when pivot is equivalent to
  // the only rotated position.
  bool is_rotation_not_needed_(const TrialSelect * select,
      const Position& pivot) {
    const Select& rotated = select->mobile();
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
