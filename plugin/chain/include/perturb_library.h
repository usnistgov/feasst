
#ifndef FEASST_CHAIN_PERTURB_LIBRARY_H_
#define FEASST_CHAIN_PERTURB_LIBRARY_H_

#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

/**
  Use an input library to orient a particle.
  In development: not fully implemented.
 */
class PerturbLibrary : public PerturbRotate {
 public:
  /**
    args:
    - library_xyz: name of xyz file that contains many configurations.
    - pivot_site: set the site index in selection with which to use as the
      pivot for rotation (default: 0).
   */
  explicit PerturbLibrary(argtype args = argtype());
  explicit PerturbLibrary(argtype * args);

  /// Read library_xyz and store the particle configuration.
  void precompute(TrialSelect * select, System * system) override;

  /// Library of particle positions.
  const std::vector<std::vector<Position> >& xyz() const { return xyz_; }

  /// Randomly select a configuration from a library,
  /// then rotate about TrialSelectParticle::site.
  void move(const bool is_position_held, System * system, TrialSelect * select,
    Random * random) override;
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbLibrary(std::istream& istr);
  virtual ~PerturbLibrary() {}

 protected:
  void serialize_perturb_library_(std::ostream& ostr) const;

 private:
  std::string library_xyz_file_name_;
  std::vector<std::vector<Position> > xyz_;
  int pivot_site_;
};

inline std::shared_ptr<PerturbLibrary> MakePerturbLibrary(
    argtype args = argtype()) {
  return std::make_shared<PerturbLibrary>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_LIBRARY_H_
