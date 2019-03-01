
#ifndef FEASST_FLAT_HISTOGRAM_BIAS_TRANSITION_MATRIX_H_
#define FEASST_FLAT_HISTOGRAM_BIAS_TRANSITION_MATRIX_H_

#include "flat_histogram/include/bias.h"

namespace feasst {

class BiasTransitionMatrix : public Bias {
 public:
  virtual ~BiasTransitionMatrix() {}
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_BIAS_TRANSITION_MATRIX_H_
