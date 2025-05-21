#include "utils/test/utils.h"

namespace feasst {

TEST(MACROS, FEASST_VERSION) {
  std::cout << FEASST_VERSION << std::endl;
}

TEST(MACROS, FEASST_INSTALL_DIR) {
  std::cout << FEASST_INSTALL_DIR << std::endl;
}

}  // namespace feasst
