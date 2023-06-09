
#ifndef FEASST_NETCDF_MOVIE_NETCDF_H_
#define FEASST_NETCDF_MOVIE_NETCDF_H_

#include "monte_carlo/include/analyze.h"
#include "netcdf/include/file_netcdf.h"

namespace feasst {

/**
  Write a trajectory of the site positions using FileXYZ format.
  Does not overwrite existing file by default.
 */
class MovieNETCDF : public AnalyzeWriteOnly {
 public:
  /**
    args:
    - append is always set to true via Stepper:set_append().
    - FileNETCDF arguments (e.g., group_index).
    - Analyze and Stepper base class arguments.
   */
  explicit MovieNETCDF(argtype args = argtype());
  explicit MovieNETCDF(argtype * args);

  /// Write the sample VMD files and the initial configuration.
  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  /// Write the configuration.
  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("MovieNETCDF"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<MovieNETCDF>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<MovieNETCDF>(args); }
  MovieNETCDF(std::istream& istr);

 private:
  std::shared_ptr<FileNETCDF> xyz_;

  std::shared_ptr<FileNETCDF> file_parse_(argtype *args) {
    args->insert({"file_name", file_name()});
    return std::make_shared<FileNETCDF>(args);
  }
};

inline std::shared_ptr<MovieNETCDF> MakeMovieNETCDF(argtype args = argtype()) {
  return std::make_shared<MovieNETCDF>(args);
}

// Placeholder class that has no arguments for depend.py to recognize this plugin.
class PlaceholderNETCDF {
 public:
  PlaceholderNETCDF() { useless_ = true; }
 private:
  bool useless_;
};

}  // namespace feasst

#endif  // FEASST_NETCDF_MOVIE_NETCDF_H_
