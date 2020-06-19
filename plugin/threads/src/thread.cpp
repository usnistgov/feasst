#include "utils/include/debug.h"
#include "threads/include/thread.h"

namespace feasst {

bool Thread::in_chunk(const int iteration, const int total_iterations) const {
  DEBUG(thread());
  DEBUG(num());
  if (iteration%num() == thread()) {
    DEBUG("true");
    return true;
  }
  DEBUG("false");
  return false;
}

}  // namespace feasst
