/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef BASE_H_
#define BASE_H_

#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <memory>
#include "./functions.h"
#include "./arguments.h"

#ifdef _OPENMP
#include <omp.h>
#endif  // _OPENMP

#ifdef MPI_H_
#include "./mpi.h"
#endif  // MPI_H_

#ifdef JSON_
#include "./../extern/json.hpp"
#endif  // JSON_

#ifdef HDF5_
#include "./H5Cpp.h"
  #ifndef H5_NO_NAMESPACE
  using namespace H5;
  #endif
#endif  // HDF5_

using std::cout;
using std::endl;
using std::vector;
using std::shared_ptr;
using std::make_shared;
using std::ofstream;
using std::ostringstream;
using std::stringstream;
using std::string;

#define STRINGIFY(FEASST_SRC_) #FEASST_SRC_
#define TOSTRING(FEASST_SRC_) STRINGIFY(FEASST_SRC_)

namespace feasst {

/**
 * Base class for all other classes to inherit.
 */
class Base {
 public:
  Base();
  virtual ~Base() {}

  /// Reconstruct object pointers after cloning (e.g., deep copy).
  void reconstruct();

  /**
   * The name of the class which should always be set at the constructor,
   * and updated for all derived classes.
   */
  std::string className() const { return className_; }

  /// The base directory where the code is located is stored in this variable.
  const char* install_dir() const { return install_dir_.c_str(); }

 protected:
  std::string className_;     //!< name of class
  string install_dir_;        //!< install directory
  int verbose_;               //!< flag for verbose printing
  Arguments argparse_;        //!< argument parser

  /// Place-holder for derived objects to perform additional reconstruction.
  virtual void reconstructDerived_() {}
};

}  // namespace feasst

#endif  // BASE_H_

