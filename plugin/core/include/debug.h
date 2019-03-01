
#ifndef FEASST_CORE_DEBUG_H_
#define FEASST_CORE_DEBUG_H_

#include <signal.h>
#include <assert.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include "core/include/custom_exception.h"

namespace feasst {

/**
  The verbosity level sets which error and debugging checks and messages
  are utilized.

  The order of verbosity is as follows:

  0. FATAL (not implemented)
  1. ERROR and ASSERT
  2. WARN and WARN_IF
  3. INFO
  4. DEBUG
  5. TRACE

  These macros will print if they are equal to or less than the verbosity level.
  Tests which expect exceptions (many) will fail if VERBOSE_LEVEL is set to 0.
  WARN_IF and ASSERT include the conditional as the first argument.

  A way to decide whether to use TRACE or DEBUG is as follows.
  If the printout is once every trial move or less frequent, then DEBUG is
  acceptable.
  If the printout is more than once every trial (e.g., inside the visitor loop
  for a model) then use TRACE instead.
*/
constexpr int VERBOSE_LEVEL = 3;

/// Used to output maximum precision to screen
/// [e.g., INFO(MAX_PRECISION << "energy: " << energy)].
#define MAX_PRECISION std::setprecision(std::numeric_limits<double>::digits10+2)

/// Return file_name with FEASST_DIR_ path removed.
std::string feasst_dir_trim_(const char* file_name);

/// If the assertion condition is not true, throw exception with message.
# define ASSERT(condition, message) \
{ \
  if (!(condition) && feasst::VERBOSE_LEVEL >= 1) { \
    std::stringstream macro_err_msg_; \
    macro_err_msg_ << "# Assertion `" #condition "` failed " \
             << feasst::feasst_dir_trim_(__FILE__) \
             << ":" << __LINE__ << ": " << message; \
    throw feasst::CustomException(macro_err_msg_); \
  } \
}

/// Throw exception with message.
# define ERROR(message) \
{ \
  if (feasst::VERBOSE_LEVEL >= 1) { \
    std::stringstream macro_err_msg_; \
    macro_err_msg_ << "# Error " \
             << feasst::feasst_dir_trim_(__FILE__) \
             << ":" << __LINE__ << ": " << message; \
    throw feasst::CustomException(macro_err_msg_); \
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

/// Warn with message to standard output if condition is true.
# define WARN_IF(condition, message) \
{ \
  if (condition && feasst::VERBOSE_LEVEL >= 2) { \
    std::stringstream macro_err_msg_; \
    std::clog << "# Warning if `" #condition "` " \
              << feasst::feasst_dir_trim_(__FILE__) \
              << ":" << __LINE__ << ": " << message << std::endl; \
  } \
}

/// Warn with message to standard output.
# define WARN(message) \
{ \
  if (feasst::VERBOSE_LEVEL >= 2) { \
    std::stringstream macro_err_msg_; \
    std::clog << "# Warning " \
              << feasst::feasst_dir_trim_(__FILE__) \
              << ":" << __LINE__ << ": " << message << std::endl; \
  } \
}

/// Inform with message to standard output.
# define INFO(message) \
{ \
  if (feasst::VERBOSE_LEVEL >= 3) { \
    std::clog << "# Info " \
              << feasst::feasst_dir_trim_(__FILE__) \
              << ":" << __LINE__ << ": " << message << std::endl; \
  } \
}

/// Debug with message to standard output.
# define DEBUG(message) \
{ \
  if (feasst::VERBOSE_LEVEL >= 4) { \
    std::clog << "# Debug " \
              << feasst::feasst_dir_trim_(__FILE__) \
              << ":" << __LINE__ << ": " << message << std::endl; \
  } \
}

/// Trace with message to standard output.
# define TRACE(message) \
{ \
  if (feasst::VERBOSE_LEVEL >= 5) { \
    std::clog << "# Trace " \
              << feasst::feasst_dir_trim_(__FILE__) \
              << ":" << __LINE__ << ": " << message << std::endl; \
  } \
}

}  // namespace feasst

#endif  // FEASST_CORE_DEBUG_H_
