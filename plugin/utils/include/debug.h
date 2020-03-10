
#ifndef FEASST_UTILS_DEBUG_H_
#define FEASST_UTILS_DEBUG_H_

#include <signal.h>
#include <assert.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>  // setprecision
#include <limits>  // numeric_limits
#ifdef _OPENMP
  #include <omp.h>
#endif  // _OPENMP
#include "utils/include/custom_exception.h"

namespace feasst {

/**
  The verbosity level sets which error and debugging checks and messages
  are utilized.

  The order of verbosity is as follows:

  0. FATAL
  1. ERROR and ASSERT
  2. WARN
  3. INFO
  4. DEBUG
  5. TRACE

  These macros will print if they are equal to or less than the verbosity level.
  Tests which expect exceptions (many) will fail if VERBOSE_LEVEL is set to 0.
  ASSERT includes the conditional as the first argument.

  A way to decide whether to use TRACE or DEBUG is as follows.
  If the printout is once every trial move or less frequent, then DEBUG is
  acceptable.
  If the printout is more than once every trial (e.g., inside the visitor loop
  for a model) then use TRACE instead.
*/
constexpr int VERBOSE_LEVEL = 3;
// HWH update verbose levels: per trial, per energy calc, per interaction
// HWH also verbosity levels for each plugin and classes.

/// Used to output maximum precision to screen
/// [e.g., INFO(MAX_PRECISION << "energy: " << energy)].
#define MAX_PRECISION std::setprecision(std::numeric_limits<double>::digits10+2)

/// Return file_name with FEASST_DIR_ path removed.
std::string feasst_dir_trim_(const char* file_name);

/// Throw exception
# ifdef _OPENMP
# define FEASST_MACRO_EXCEPTION(message, name) \
{ \
  std::stringstream macro_err_msg_; \
  macro_err_msg_ << "# " << name << omp_get_thread_num() << " " \
           << feasst::feasst_dir_trim_(__FILE__) \
           << ":" << __LINE__ << ": " << message; \
  throw feasst::CustomException(macro_err_msg_); \
}
# else  // _OPENMP
# define FEASST_MACRO_EXCEPTION(message, name) \
{ \
  std::stringstream macro_err_msg_; \
  macro_err_msg_ << "# " << name << " " \
           << feasst::feasst_dir_trim_(__FILE__) \
           << ":" << __LINE__ << ": " << message; \
  throw feasst::CustomException(macro_err_msg_); \
}
# endif  // _OPENMP

/// If the assertion condition is not true, throw exception with message.
# define ASSERT(condition, message) \
{ \
  if (!(condition) and feasst::VERBOSE_LEVEL >= 1) { \
    FEASST_MACRO_EXCEPTION(message, "Assertion `" #condition "` failed") \
  } \
}

/// Throw exception with message.
# define ERROR(message) \
{ \
  if (feasst::VERBOSE_LEVEL >= 1) { \
    FEASST_MACRO_EXCEPTION(message, "Error") \
  } \
}

/// Throw exception with message.
# define FATAL(message) \
{ \
  if (feasst::VERBOSE_LEVEL >= 0) { \
    FEASST_MACRO_EXCEPTION(message, "Fatal error") \
  } \
}

/// Try block to expect an error. Disable if VERBOSE_LEVEL is 0.
# define TRY(code) \
{ \
  if (feasst::VERBOSE_LEVEL > 0) { \
    try { \
      code \
    } \
  } \
}

/// Expect to catch exception with phrase in message.
# define CATCH_PHRASE(phrase) \
FAIL() << "Expected failure"; \
} catch(feasst::CustomException e) { \
  std::string what(e.what()); \
  std::stringstream ermsg; \
  ermsg << "The phrase(" << phrase << ") was not contained in the expected " \
    << "error message(" << what << ")"; \
  ASSERT(what.find(phrase) != std::string::npos, ermsg.str()); \
} catch(...) { \
  FAIL() << "Unrecognized exception"; \

/// Debug with message to standard output.
# ifdef _OPENMP
# define FEASST_MACRO_OUTPUT(message, name, level) \
{ \
  if (feasst::VERBOSE_LEVEL >= level) { \
    std::cout << "# " << name << omp_get_thread_num() \
              << " [" << feasst::feasst_dir_trim_(__FILE__) << ":" << __LINE__ \
              << "] " << message << std::endl; \
  } \
}
# else  // _OPENMP
# define FEASST_MACRO_OUTPUT(message, name, level) \
{ \
  if (feasst::VERBOSE_LEVEL >= level) { \
    std::cout << "# " << name \
              << " [" << feasst::feasst_dir_trim_(__FILE__) << ":" << __LINE__ \
              << "] " << message << std::endl; \
  } \
}
#endif  // _OPENMP

/// Define and name the various levels of verbosity.
# define WARN(message) FEASST_MACRO_OUTPUT(message, "Warn ", 2)
# define INFO(message) FEASST_MACRO_OUTPUT(message, "Info ", 3)
# define DEBUG(message) FEASST_MACRO_OUTPUT(message, "Debug", 4)
# define TRACE(message) FEASST_MACRO_OUTPUT(message, "Trace", 5)
}  // namespace feasst

#endif  // FEASST_UTILS_DEBUG_H_
