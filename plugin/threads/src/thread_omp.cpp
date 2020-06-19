#ifdef _OPENMP
  #include <omp.h>
#endif // _OPENMP
#include "utils/include/debug.h"
#include "threads/include/thread_omp.h"

namespace feasst {

ThreadOMP::ThreadOMP() {
  num_ = 1;
  thread_ = 0;
  if (!is_enabled()) WARN("OMP is not enabled");
  #ifdef _OPENMP
    num_ = static_cast<int>(omp_get_num_threads());
    thread_ = omp_get_thread_num();
  #endif // _OPENMP
}

bool ThreadOMP::is_enabled() const {
#ifdef _OPENMP
  return true;
#else // _OPENMP
  return false;
#endif // _OPENMP
}

}  // namespace feasst
