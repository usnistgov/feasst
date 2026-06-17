
#ifndef FEASST_ACTION_ENERGY_GRID_H_
#define FEASST_ACTION_ENERGY_GRID_H_

#include <memory>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Generate a grid of energies for visualization, plotting and testing.
 */
class EnergyGrid : public Action {
 public:
  //@{
  /** @name Arguments
    - config: Name of Configuration (default: 0).
    - output_file: Name of file to output energy grid.
    - resolution: Comma-separated values of ranges for each dimension.
      Ranges are defined as colon-separated min:max:increment values.
   */
  explicit EnergyGrid(argtype args = argtype());
  explicit EnergyGrid(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<EnergyGrid>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<EnergyGrid>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit EnergyGrid(std::istream& istr);
  virtual ~EnergyGrid();

  //@}
 private:
  std::string config_, output_file_;
  std::vector<std::vector<double> > res_;
};

}  // namespace feasst

#endif  // FEASST_ACTION_ENERGY_GRID_H_
