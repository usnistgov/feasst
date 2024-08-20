
#ifndef FEASST_MPI_MODEL_MPI_H_
#define FEASST_MPI_MODEL_MPI_H_

#include <string>
#include <memory>
#include "system/include/model_two_body.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class ThreadMPI;

/**
  Establish an MPI interface to obtain the interaction potential from another
  thread.
  The other thread receives the following as bytes which are the following three

  squared_distance, type1, type2

  where the squared_distance is the square of the distance between the centers
  of two sites, and the types are the types of the sites (see Configuration).

  The MPI client which returns the energy receives a negative squared distance
  when it is no longer needed.
  This happens upon the destruction of ModelMPI.
 */
class ModelMPI : public ModelTwoBody {
 public:
  //@{
  /** @name Arguments
    - hard_sphere_threshold: when r < threshold*sigma, return NEAR_INFINITY
      (default: 0.2).
   */
  explicit ModelMPI(argtype args = argtype());
  explicit ModelMPI(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  //void precompute(const ModelParams& existing) override;
  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  /// Return the threshold for hard sphere interaction.
  double hard_sphere_threshold() const;
  const double& hard_sphere_threshold_sq() const {
    return hard_sphere_threshold_sq_; }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelMPI>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ModelMPI>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ModelMPI(std::istream& istr);
  virtual ~ModelMPI();

  //@}
 protected:
  void serialize_model_mpi_(std::ostream& ostr) const;

 private:
  double hard_sphere_threshold_sq_;

  // initialization and serialization not required
  double en_;
  std::shared_ptr<ThreadMPI> thread_;
};

inline std::shared_ptr<ModelMPI> MakeModelMPI(
    argtype args = argtype()) {
  return std::make_shared<ModelMPI>(args);
}

}  // namespace feasst

#endif  // FEASST_MPI_MODEL_MPI_H_
