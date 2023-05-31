#include "utils/test/utils.h"
#include "fftw/include/scattering_fftw.h"

namespace feasst {

TEST(ScatteringFFTW, serialize) {
  auto an = MakeScatteringFFTW();
  auto an2 = test_serialize<ScatteringFFTW, Analyze>(*an);
}

}  // namespace feasst
