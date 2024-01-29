#include "mpi.h"
#include "utils/include/debug.h"
#include "mpi/include/thread_mpi.h"

namespace feasst {

ThreadMPI::ThreadMPI() {
  num_ = 1;
  thread_ = 0;
  int initialized;
  MPI_Initialized(&initialized);
  if (!initialized) {
    MPI_Init(NULL, NULL);
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &thread_);
  MPI_Comm_size(MPI_COMM_WORLD, &num_);
  DEBUG("thread: " << thread_);
  DEBUG("num: " << num_);
}

ThreadMPI::~ThreadMPI() {
  int finalized;
  MPI_Finalized(&finalized);
  if (!finalized) {
     MPI_Finalize();
  }
}

bool ThreadMPI::is_enabled() const {
  return true;
}

}  // namespace feasst
