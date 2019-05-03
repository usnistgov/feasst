
#ifndef FEASST_CORE_MOVIE_H_
#define FEASST_CORE_MOVIE_H_

#include "core/include/analyze.h"

namespace feasst {

/**
  Write a trajectory of the site positions using XYZ file format.
  Does not overwrite existing file.
 */
// HWH cite XYZ.
class Movie : public AnalyzeWriteOnly {
 public:
  Movie(const argtype &args = argtype()) : AnalyzeWriteOnly(args) {
    set_append();

    // require file_name
    args_.init(args);
    ASSERT(args_.key("file_name").used(), "file name is required");
  }

  /// Write the sample VMD files and the initial configuration.
  void initialize(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    ASSERT(!file_name().empty(), "file name required. Did you forget to " <<
      "Analyze::set_file_name()?");

    // write xyz
    xyz_.set_append(1);
    if (state() == criteria->state()) {
      xyz_.write(file_name(), system.configuration());
    }

    // write vmd
    std::stringstream ss;
    ss << file_name() << ".vmd";
    vmd_.write(ss.str(), system.configuration(), file_name());
  }

  /// Write the configuration.
  std::string write(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    // ensure the following order matches the header from initialization.
    xyz_.write(file_name(), system.configuration());
    return std::string("");
  }

  std::string class_name() const override { return std::string("Movie"); }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
    feasst_serialize_fstobj(xyz_, ostr);
    feasst_serialize_fstobj(vmd_, ostr);
  }

  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Movie>(istr); }

  Movie(std::istream& istr) : AnalyzeWriteOnly(istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize_fstobj(&xyz_, istr);
    feasst_deserialize_fstobj(&vmd_, istr);
  }

 private:
  FileXYZ xyz_;
  FileVMD vmd_;
};

inline std::shared_ptr<Movie> MakeMovie(const argtype &args = argtype()) {
  return std::make_shared<Movie>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_MOVIE_H_
