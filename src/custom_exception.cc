/**
 * \file
 *
 * \brief randomly selects monte carlo trials
 *
 */

#ifdef MPI_H_
  #include <mpi.h>
#endif
#ifdef _OPENMP
  #include <omp.h>
#endif  // _OPENMP
#include "./custom_exception.h"

customException::customException(std::string m
  ) {
  msg_ = m;
  catMessage();
  cout << msg_ << endl;
  exit(0);  // terminate immediately for easy backtrace
}
customException::customException(std::stringstream& m) {
  msg_ = m.str();
  catMessage();
  cout << msg_ << endl;
  exit(0);  // terminate immediately for easy backtrace
}

void customException::catMessage() {
  int nproc = 0;
  #ifdef MPI_H_
    int initialized;
    MPI_Initialized(&initialized);
    if (initialized) MPI_Comm_rank(MPI_COMM_WORLD, &nproc);
  #endif  // MPI_H_
  #ifdef _OPENMP
    nproc = omp_get_thread_num();
    #pragma omp critical
    {
  #endif  // _OPENMP

  std::stringstream message;
  message << "Throw on proc " << nproc << " : " << msg_;
  msg_.assign(message.str());

  #ifdef _OPENMP
    }
  #endif  // _OPENMP

}
