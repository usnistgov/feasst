
#ifndef FEASST_THREADS_THREADS_OMP_H_
#define FEASST_THREADS_THREADS_OMP_H_

#include <memory>
#include "threads/include/threads.h"

namespace feasst {

class ThreadsOMP : public Threads {
 public:
  bool is_enabled() const override;
  int num() const override;
  virtual ~ThreadsOMP() {}
};

inline std::shared_ptr<ThreadsOMP> MakeThreadsOMP() {
  return std::make_shared<ThreadsOMP>();
}

}  // namespace feasst

#endif  // FEASST_THREADS_THREADS_OMP_H_
