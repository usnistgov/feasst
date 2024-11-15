
#ifndef FEASST_ANISO_BACKMAP_H_
#define FEASST_ANISO_BACKMAP_H_

#include "configuration/include/file_vmd.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Write a backmapping trajector where anisotropic sites are converted into
  multi-site models.
 */
class Backmap : public AnalyzeWriteOnly {
 public:
  //@{
  /** @name Arguments
    - site[i]: site type of anisotropic site to backmap.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
    - fstprt[i]: fstprt to backmap anisotropic site into.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
    - Stepper arguments.
    - append is always set to true via Stepper:set_append().
   */
  explicit Backmap(argtype args = argtype());
  explicit Backmap(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Write the sample VMD files and the initial configuration.
  void initialize(MonteCarlo * mc) override;

  /// Write the configuration.
  std::string write(const MonteCarlo& mc) override;

  // serialize
  std::string class_name() const override { return std::string("Backmap"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Backmap>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<Backmap>(args); }
  Backmap(std::istream& istr);

  //@}
 private:
  std::vector<int> site_types_;
  std::vector<std::string> site_fstprt_;
  std::shared_ptr<Configuration> all_atom_;
  FileXYZ xyz_;
  FileVMD vmd_;

  void add_backmap_particles_();
};

}  // namespace feasst

#endif  // FEASST_ANISO_BACKMAP_H_
