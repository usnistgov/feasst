#include "utils/test/utils.h"
#include "utils/include/utils_io.h"
#include "threads/include/thread_omp.h"

namespace feasst {

TEST(ThreadOMP, chunk) {
  const int num = 33;
  std::vector<int> data(num, -1);
  #ifdef _OPENMP
  #pragma omp parallel
  {
    auto thread = MakeThreadOMP();
    for (int i = 0; i < num; ++i) {
      if (thread->in_chunk(i, num)) {
        EXPECT_EQ(data[i], -1);
        data[i] = thread->thread();
      }
    }
  }
  #endif // _OPENMP

  INFO(feasst_str(data));
}

}  // namespace feasst
