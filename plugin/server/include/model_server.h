
#ifndef FEASST_SERVER_MODEL_SERVER_H_
#define FEASST_SERVER_MODEL_SERVER_H_

#include <string>
#include <memory>
#include "system/include/model_two_body.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class Server;

/**
  Establish a server to obtain the interaction potential from a client.
  The client receives byte strings which are the following three csv,

  squared_distance, type1, type2

  where the squared_distance is the square of the distance between the centers
  of two sites, and the types are the types of the sites (see Configuration).
 */
class ModelServer : public ModelTwoBody {
 public:
  //@{
  /** @name Arguments
    - hard_sphere_threshold: when r < threshold*sigma, return NEAR_INFINITY
      (default: 0.2).
    - Server arguments.
   */
  explicit ModelServer(argtype args = argtype());
  explicit ModelServer(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

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
    return std::make_shared<ModelServer>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ModelServer>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ModelServer(std::istream& istr);
  virtual ~ModelServer();

  //@}
 protected:
  void serialize_model_server_(std::ostream& ostr) const;

 private:
  double hard_sphere_threshold_sq_;
  std::unique_ptr<Server> server_;
};

}  // namespace feasst

#endif  // FEASST_SERVER_MODEL_SERVER_H_
