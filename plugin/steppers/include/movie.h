
#ifndef FEASST_STEPPERS_MOVIE_H_
#define FEASST_STEPPERS_MOVIE_H_

#include "configuration/include/file_vmd.h"
#include "configuration/include/file_xyz.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

// HWH allow different formats.
// HWH for example, incorportate FileXYZPatch
/**
  Write a trajectory of the site positions using FileXYZ format.
  Appends to existing file by default.
 */
class Movie : public AnalyzeWriteOnly {
 public:
  //@{
  /** @name Arguments
    - FileXYZ arguments (e.g., group_index).
    - FileVMD arguments (e.g., min_sigma).
    - Stepper arguments.
    - append is always set to true via Stepper:set_append().
   */
  explicit Movie(argtype args = argtype());
  explicit Movie(argtype * args);

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
  std::string class_name() const override { return std::string("Movie"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Movie>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<Movie>(args); }
  Movie(std::istream& istr);

  //@}
 private:
  FileXYZ xyz_;
  FileVMD vmd_;
};

inline std::shared_ptr<Movie> MakeMovie(argtype args = argtype()) {
  return std::make_shared<Movie>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_MOVIE_H_
