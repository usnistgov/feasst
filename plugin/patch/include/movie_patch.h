
#ifndef FEASST_PATCH_MOVIE_PATCH_H_
#define FEASST_PATCH_MOVIE_PATCH_H_

#include "monte_carlo/include/analyze.h"
#include "patch/include/file_xyz_patch.h"

namespace feasst {

/**
  Write a trajectory of the site positions using FileXYZPatch format.
  Appends to existing file by default.
 */
class MoviePatch : public AnalyzeWriteOnly {
 public:
  //@{
  /** @name Arguments
    - FileXYZPatch arguments.
    - Stepper arguments.
    - append is always set to true via Stepper:set_append().
   */
  explicit MoviePatch(argtype args = argtype());
  explicit MoviePatch(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Write the sample VMD files and the initial configuration.
  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  /// Write the configuration.
  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("MoviePatch"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<MoviePatch>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<MoviePatch>(args); }
  MoviePatch(std::istream& istr);

  //@}
 private:
  FileXYZPatch xyz_;
  FileVMDPatch vmd_;
};

inline std::shared_ptr<MoviePatch> MakeMoviePatch(argtype args = argtype()) {
  return std::make_shared<MoviePatch>(args);
}

}  // namespace feasst

#endif  // FEASST_PATCH_MOVIE_PATCH_H_
