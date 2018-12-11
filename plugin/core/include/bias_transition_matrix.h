
#ifndef FEASST_CORE_BIAS_TRANSITION_MATRIX_H_
#define FEASST_CORE_BIAS_TRANSITION_MATRIX_H_

#include "core/include/bias.h"

namespace feasst {

class BiasTransitionMatrix : public Bias {
 public:
  virtual ~BiasTransitionMatrix() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_BIAS_TRANSITION_MATRIX_H_
