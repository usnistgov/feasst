
#ifndef FEASST_THREAD_THREAD_OMP_H_
#define FEASST_THREAD_THREAD_OMP_H_

#include <memory>
#include "threads/include/thread.h"

namespace feasst {

class ThreadOMP : public Thread {
 public:
  ThreadOMP();
  bool is_enabled() const override;
  int num() const override { return num_; }
  int thread() const override { return thread_; }
  virtual ~ThreadOMP() {}

 private:
  int num_;
  int thread_;
};

inline std::shared_ptr<ThreadOMP> MakeThreadOMP() {
  return std::make_shared<ThreadOMP>();
}

}  // namespace feasst

#endif  // FEASST_THREAD_THREAD_OMP_H_
