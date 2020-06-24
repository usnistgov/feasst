#include "utils/include/debug.h"
#include "threads/include/thread.h"

namespace feasst {

bool Thread::in_chunk(const int iteration,
    const int total_iterations,
    const int node,
    const int num_nodes) const {
  DEBUG(thread());
  DEBUG(num());
  if (iteration%(num_nodes*num()) == thread() + node*num()) {
    DEBUG("true");
    return true;
  }
  DEBUG("false");
  return false;
}

}  // namespace feasst
