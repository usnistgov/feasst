#include "utils/test/utils.h"
#include "utils/include/io.h"
#include "mpi/include/thread_mpi.h"

namespace feasst {

TEST(ThreadMPI, chunk_LONG) {
  const int num = 33;
  std::vector<int> data(num, -1);
  auto thread = MakeThreadMPI();
  for (int i = 0; i < num; ++i) {
    if (thread->in_chunk(i, num)) {
      EXPECT_EQ(data[i], -1);
      data[i] = thread->thread();
    }
  }
}

}  // namespace feasst
