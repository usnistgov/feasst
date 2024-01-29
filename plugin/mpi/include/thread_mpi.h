
#ifndef FEASST_MPI_THREAD_MPI_H_
#define FEASST_MPI_THREAD_MPI_H_

#include "mpi.h"
#include <memory>
#include "threads/include/thread.h"

namespace feasst {

class ThreadMPI : public Thread {
 public:
  ThreadMPI();
  bool is_enabled() const override;
  int num() const override { return num_; }
  int thread() const override { return thread_; }
  virtual ~ThreadMPI();

 private:
  int num_;
  int thread_;
};

inline std::shared_ptr<ThreadMPI> MakeThreadMPI() {
  return std::make_shared<ThreadMPI>();
}

}  // namespace feasst

#endif  // FEASST_MPI_THREAD_MPI_H_
