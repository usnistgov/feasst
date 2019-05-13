#include <gtest/gtest.h>
#include "utils/include/debug.h"

namespace feasst {

TEST(Debug, ASSERT) {
  try {
    ASSERT(0, "failure");
    CATCH_PHRASE("failure");
  }
}

TEST(Debug, ERROR) {
  try {
    ERROR("failure");
    CATCH_PHRASE("failure");
  }
}

TEST(Debug, WARN_INFO_DEBUG_TRACE) {
  WARN_IF(1 > 0, "WARN_IF");
  WARN("WARN");
  INFO("INFO");
  DEBUG("DEBUG");
  TRACE("TRACE");
}

TEST(Debug, VERBOSE_LEVEL) {
  const int threshold = 3;
  if (VERBOSE_LEVEL < threshold) {
    std::cout << "By default, we expect a verbosity level of atleast "
              << threshold << ". This test failure serves as a warning."
              << std::endl;
  }
  EXPECT_GE(VERBOSE_LEVEL, threshold);
}

TEST(Debug, feasst_dir_trim_) {
  std::string dir = feasst_dir_trim_(FEASST_DIR_);
  EXPECT_EQ(0, dir.size());
}

}  // namespace feasst
