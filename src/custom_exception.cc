/**
 * \file
 *
 * \brief randomly selects monte carlo trials
 *
 */

#ifdef MPI_H_
  #include <mpi.h>
#endif
#ifdef OMP_H_
  #include <omp.h>
#endif  // OMP_H_
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
}

customException::~customException() {
  throw(msg_);
}

void customException::catMessage() {
  int nproc = 0;
  #ifdef MPI_H_
    int initialized;
    MPI_Initialized(&initialized);
    if (initialized) MPI_Comm_rank(MPI_COMM_WORLD, &nproc);
  #endif  // MPI_H_
  #ifdef OMP_H_
    nproc = omp_get_thread_num();
    #pragma omp critical
    {
  #endif  // OMP_H_

  std::stringstream message;
  message << "Throw on proc " << nproc << " : " << msg_;
  msg_.assign(message.str());

  #ifdef OMP_H_
    }
  #endif  // OMP_H_

}
