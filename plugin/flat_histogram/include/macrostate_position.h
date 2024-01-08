
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_POSITION_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_POSITION_H_

#include "flat_histogram/include/macrostate.h"

namespace feasst {

/**
  Defines the macrostate to be the position of a site on a particule in a
  dimension.
 */
class MacrostatePosition : public Macrostate {
 public:
  //@{
  /** @name Arguments
    - particle_index: index of particle (default: 0).
    - site_index: index of site (default: 0).
    - dimension: dimension of position (default: 0).
    - Macrostate arguments.
  */
  explicit MacrostatePosition(argtype args = argtype());
  explicit MacrostatePosition(argtype * args) :
    MacrostatePosition(Histogram(args), args) {}

  //@}
  /** @name Public Functions
   */
  //@{

  /// Arguments as described above, but with explicit histogram object.
  MacrostatePosition(const Histogram& histogram, argtype args = argtype());
  MacrostatePosition(const Histogram& histogram, argtype * args);

  double value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override;
  std::shared_ptr<Macrostate> create(std::istream& istr) const override;
  std::shared_ptr<Macrostate> create(argtype * args) const override;
  void serialize(std::ostream& ostr) const override;
  MacrostatePosition(std::istream& istr);
  virtual ~MacrostatePosition() {}

  //@}
 private:
  int particle_index_;
  int site_index_;
  int dimension_;
};

inline std::shared_ptr<MacrostatePosition> MakeMacrostatePosition(
    const Histogram& histogram, argtype args = argtype()) {
  return std::make_shared<MacrostatePosition>(histogram, args);
}

inline std::shared_ptr<MacrostatePosition> MakeMacrostatePosition(
    argtype args = argtype()) {
  return std::make_shared<MacrostatePosition>(args);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_MACROSTATE_POSITION_H_
