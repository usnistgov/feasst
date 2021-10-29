#include "threads/include/thread.h"

namespace feasst {

bool Thread::in_chunk(const int iteration,
    const int total_iterations,
    const int node,
    const int num_nodes) const {
  if (iteration%(num_nodes*num()) == thread() + node*num()) {
    return true;
  }
  return false;
}

}  // namespace feasst
