/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include <iostream>
#ifdef MPI_H_
  #include <mpi.h>
#endif
#ifdef _OPENMP
  #include <omp.h>
#endif  // _OPENMP
#include "utils/include/custom_exception.h"

namespace feasst {

CustomException::CustomException(std::stringstream& m) {
  msg_ = m.str();
  catMessage_();

  // Output error message upon construction if openMP is used.
  // Otherwise, you may never see the error message!
  #ifdef _OPENMP
    if (omp_get_num_threads() > 1) {
      //#pragma omp critical
     // {
        std::cout << msg_ << std::endl;
     // }
    }
  #endif  // _OPENMP
}

void CustomException::catMessage_() {
  int nproc = 0;
  #ifdef MPI_H_
    int initialized;
    MPI_Initialized(&initialized);
    if (initialized) MPI_Comm_rank(MPI_COMM_WORLD, &nproc);
  #endif  // MPI_H_
  #ifdef _OPENMP
    nproc = omp_get_thread_num();
//    #pragma omp critical
//    {
  #endif  // _OPENMP

  std::stringstream message;
  message << "Throw on proc " << nproc << " : " << msg_;
  msg_.assign(message.str());

//  #ifdef _OPENMP
//    }
//  #endif  // _OPENMP
}

}  // namespace feasst
