
#ifndef FEASST_CORE_MOVIE_H_
#define FEASST_CORE_MOVIE_H_

#include "core/include/analyze.h"

namespace feasst {

class Movie : public AnalyzeWriteOnly {
 public:
  Movie(const argtype &args = argtype()) : AnalyzeWriteOnly(args) {
    set_append();
  }
  void initialize(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    ASSERT(!file_name().empty(), "file name required. Did you forget to " <<
      "Analyze::set_file_name()?");
    xyz_.write(file_name(), system.configuration());
    xyz_.set_append(1);

    // write vmd
    std::stringstream ss;
    ss << file_name() << ".vmd";
    vmd_.write(ss.str(), system.configuration(), file_name());
  }

  std::string write(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    // ensure the following order matches the header from initialization.
    xyz_.write(file_name(), system.configuration());
    return std::string("");
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
    feasst_serialize_fstobj(xyz_, ostr);
    feasst_serialize_fstobj(vmd_, ostr);
  }

  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    feasst_deserialize_version(istr);
    auto analyze = std::make_shared<Movie>();
    feasst_deserialize_fstobj(&(analyze->xyz_), istr);
    feasst_deserialize_fstobj(&(analyze->vmd_), istr);
    return analyze;
  }

 private:
  const std::string class_name_ = "Movie";
  FileXYZ xyz_;
  FileVMD vmd_;
};

inline std::shared_ptr<Movie> MakeMovie(const argtype &args = argtype()) {
  return std::make_shared<Movie>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_MOVIE_H_
