#include "utils/test/utils.h"
#include "aniso/include/combine_parallel_files.h"

namespace feasst {

TEST(CombineParallelFiles, serialize) {
  CombineParallelFiles obj;
  CombineParallelFiles obj2 = test_serialize(obj);
}

}  // namespace feasst
