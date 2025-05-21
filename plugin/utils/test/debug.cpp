#include <sstream>
#include "utils/test/utils.h"
#include "utils/include/debug.h"

namespace feasst {

TEST(Debug, ASSERT) {
  TRY(
    ASSERT(0, "failure");
    CATCH_PHRASE("failure");
  );
}

TEST(Debug, ERROR) {
  TRY(
    ERROR("failure");
    CATCH_PHRASE("failure");
  );
}

TEST(Debug, WARN_INFO_DEBUG_TRACE) {
  WARN("WARN");
  INFO("INFO");
  DEBUG("DEBUG");
  TRACE("TRACE");
}

TEST(Debug, feasst_dir_trim_) {
  std::stringstream ss;
  ss << FEASST_INSTALL_DIR << "/";
  std::string dir = feasst_dir_trim_(ss.str().c_str());
  std::cout << "# dir:" << dir << std::endl;
  EXPECT_EQ(0, dir.size());
}

}  // namespace feasst
