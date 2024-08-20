
#ifndef FEASST_SERVER_VISIT_MODEL_INNER_SERVER_H_
#define FEASST_SERVER_VISIT_MODEL_INNER_SERVER_H_

#include <memory>
#include "math/include/matrix.h"
#include "math/include/euler.h"
#include "system/include/visit_model_inner.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class Server;

/**
  Model anisotropic sites by the relative orientation of five angles in 3D
  as described in VisitModelInnerTable.

  The client receives byte strings which are the following eight csv,

  squared_distance, s1, s2, e1, e2, e3, type1, type2

  where the squared_distance is the square of the distance between the centers
  of two sites, the types are the types of the sites (see Configuration), the
  s1, s2 are the spherical coordinate angles as descirbed in the Position class,
  and e1, e2, e3 are the Euler angles.
 */
class VisitModelInnerServer : public VisitModelInner {
 public:
  //@{
  /** @name Arguments
    - ignore_energy: do not read the energy table (default: false).
    - server_site[i]: add the i-th type of site to include with this potential.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
    - Server arguments.
   */
  explicit VisitModelInnerServer(argtype args = argtype());
  explicit VisitModelInnerServer(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(Configuration * config) override;
  void compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    ModelTwoBody * model,
    const bool is_old_config,
    Position * relative,
    Position * pbc,
    const double weight = 1.) override;

  std::shared_ptr<VisitModelInner> create(std::istream& istr) const override {
    return std::make_shared<VisitModelInnerServer>(istr); }
  std::shared_ptr<VisitModelInner> create(argtype * args) const override {
    return std::make_shared<VisitModelInnerServer>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelInnerServer(std::istream& istr);
  virtual ~VisitModelInnerServer() {}

  //@}
 private:
  int aniso_index_ = -1;
  std::vector<int> site_types_;
  std::vector<int> t2index_;
  std::unique_ptr<Server> server_;

  // no serialized optimization variables
  Position pos1_, pos2_, sph_;
  RotationMatrix rot1_, rot2_, rot3_;
  Euler euler_;
};

inline std::shared_ptr<VisitModelInnerServer> MakeVisitModelInnerServer(
    argtype args = argtype()) {
  return std::make_shared<VisitModelInnerServer>(args);
}

}  // namespace feasst

#endif  // FEASST_SERVER_VISIT_MODEL_INNER_SERVER_H_
