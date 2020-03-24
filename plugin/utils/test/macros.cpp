#include "utils/test/utils.h"

namespace feasst {

TEST(MACROS, FEASST_VERSION) {
  std::cout << FEASST_VERSION_ << std::endl;
}

TEST(MACROS, FEASST_DIR) {
  std::cout << FEASST_DIR_ << std::endl;
}

}  // namespace feasst
