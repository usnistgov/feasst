
#ifndef FEASST_STEPPERS_MOVIE_H_
#define FEASST_STEPPERS_MOVIE_H_

#include "configuration/include/file_xyz.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

// HWH cite XYZ or allow different formats.
/**
  Write a trajectory of the site positions using XYZ file format.
  Does not overwrite existing file by default.
 */
class Movie : public AnalyzeWriteOnly {
 public:
  Movie(const argtype &args = argtype());

  /// Write the sample VMD files and the initial configuration.
  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  /// Write the configuration.
  std::string write(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("Movie"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Movie>(istr); }
  Movie(std::istream& istr);

 private:
  FileXYZ xyz_;
  FileVMD vmd_;
};

inline std::shared_ptr<Movie> MakeMovie(const argtype &args = argtype()) {
  return std::make_shared<Movie>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_MOVIE_H_
