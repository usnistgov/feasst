#include "utils/test/utils.h"
#include "models/include/dihedral_ryckaert_bellemans.h"

namespace feasst {

TEST(DihedralRyckaertBellemans, serialize) {
  DihedralRyckaertBellemans dihedral;
  auto dihedral2 = test_serialize(dihedral);
}

}  // namespace feasst
