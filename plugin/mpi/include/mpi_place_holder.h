
#ifndef FEASST_MPI_MPI_PLACE_HOLDER_H_
#define FEASST_MPI_MPI_PLACE_HOLDER_H_

#include <string>
#include <memory>
#include "utils/include/arguments.h"
#include "system/include/model_two_body.h"

namespace feasst {

/**
 */
class MPIPlaceHolder : public ModelTwoBody {
 public:
  explicit MPIPlaceHolder(argtype args = argtype());
  explicit MPIPlaceHolder(argtype * args);

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<MPIPlaceHolder>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<MPIPlaceHolder>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit MPIPlaceHolder(std::istream& istr);
  virtual ~MPIPlaceHolder() {}
};

inline std::shared_ptr<MPIPlaceHolder> MakeMPIPlaceHolder(
    argtype args = argtype()) {
  return std::make_shared<MPIPlaceHolder>(args);
}

}  // namespace feasst

#endif  // FEASST_MPI_MPI_PLACE_HOLDER_H_
