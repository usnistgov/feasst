
#ifndef FEASST_SPHEROCYLINDER_MOVIE_SPHEROCYLINDER_H_
#define FEASST_SPHEROCYLINDER_MOVIE_SPHEROCYLINDER_H_

#include "monte_carlo/include/analyze.h"
#include "patch/include/file_xyz_spherocylinder.h"

namespace feasst {

/**
  Write a trajectory of the site positions using FileXYZSpherocylinder format.
  Appends to existing file by default.
 */
class MovieSpherocylinder : public AnalyzeWriteOnly {
 public:
  //@{
  /** @name Arguments
    - FileXYZSpherocylinder arguments.
    - Stepper arguments.
    - append is always set to true via Stepper:set_append().
   */
  explicit MovieSpherocylinder(argtype args = argtype());
  explicit MovieSpherocylinder(argtype * args);

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
  std::string class_name() const override { return std::string("MovieSpherocylinder"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<MovieSpherocylinder>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<MovieSpherocylinder>(args); }
  MovieSpherocylinder(std::istream& istr);

  //@}
 private:
  FileXYZSpherocylinder xyz_;
  FileVMDSpherocylinder vmd_;
};

inline std::shared_ptr<MovieSpherocylinder> MakeMovieSpherocylinder(argtype args = argtype()) {
  return std::make_shared<MovieSpherocylinder>(args);
}

}  // namespace feasst

#endif  // FEASST_SPHEROCYLINDER_MOVIE_SPHEROCYLINDER_H_
