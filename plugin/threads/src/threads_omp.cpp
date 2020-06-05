#ifdef _OPENMP
  #include <omp.h>
#endif // _OPENMP
#include "utils/include/debug.h"
#include "threads/include/threads_omp.h"

namespace feasst {

bool ThreadsOMP::is_enabled() const {
#ifdef _OPENMP
  return true;
#else // _OPENMP
  return false;
#endif // _OPENMP
}

int ThreadsOMP::num() const {
  int num_threads = 1;
  ASSERT(is_enabled(), "OMP is not enabled");
  #ifdef _OPENMP
    #pragma omp parallel
    {
      int proc_id = omp_get_thread_num();
      DEBUG("hello from " << proc_id);
      if (proc_id == 0) {
        num_threads = static_cast<int>(omp_get_num_threads());
      }
    }
  #endif // _OPENMP
  return num_threads;
}

}  // namespace feasst
