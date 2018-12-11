
#ifndef FEASST_EXAMPLE_MODEL_EXAMPLE_H_
#define FEASST_EXAMPLE_MODEL_EXAMPLE_H_

#include "core/include/model_two_body.h"

namespace feasst {

/**
  Documentation for new class.
 */
class ModelExample : public ModelTwoBody {
 public:
  // Return the energy between sites with a given relative position.
  double evaluate(const Position &relative,
                  const Site& site1,
                  const Site& site2,
                  const ModelParams& model_params) const;

  virtual ~ModelExample() {}

 private:
};

}  // namespace feasst

#endif  // FEASST_EXAMPLE_MODEL_EXAMPLE_H_
