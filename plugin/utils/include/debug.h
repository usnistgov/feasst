
#ifndef FEASST_UTILS_DEBUG_H_
#define FEASST_UTILS_DEBUG_H_

#include <string>
#include <sstream>
#include "utils/include/definitions.h"
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

  A way to decide whether to use TRACE or DEBUG is as follows.
  If the printout is once every trial move or less frequent, then DEBUG is
  acceptable.
  If the printout is more than once every trial (e.g., inside the visitor loop
  for a model) then use TRACE instead.

  Set the VERBOSE_LEVEL during installation with:
  "cmake -DFEASST_VERBOSE_LEVEL=3 .."

  The VERBOSE_LEVEL is set in utils/include/defintions.h
*/
// HWH update verbose levels: per trial, per energy calc, per interaction
// HWH also verbosity levels for each plugin and classes.

/// Return file_name with FEASST_DIR_ path removed.
std::string feasst_dir_trim_(const char* file_name);

/// Throw exception
# define FEASST_MACRO_EXCEPTION(message, name) \
{ \
  std::stringstream macro_err_msg_; \
  macro_err_msg_ << "#" << name << feasst::feasst_omp_thread() << " " \
           << feasst::feasst_dir_trim_(__FILE__) \
           << ":" << __LINE__ << ": " << message; \
  throw feasst::CustomException(macro_err_msg_); \
}

std::string feasst_omp_thread();

/// If the assertion condition is not true, throw exception with message.
# define ASSERT(condition, message) \
{ \
  if (!(condition) && feasst::VERBOSE_LEVEL >= 1) { \
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
FAIL() << "Did not find expected failure."; \
} catch(const feasst::CustomException& e) { \
  std::string what(e.what()); \
  std::stringstream ermsg; \
  ermsg << "The phrase(" << phrase << ") was not contained in the expected " \
    << "error message(" << what << ")"; \
  ASSERT(what.find(phrase) != std::string::npos, ermsg.str()); \
} catch(...) { \
  FAIL() << "Unrecognized exception"; \

/// Debug with message to standard output.
# define FEASST_MACRO_OUTPUT(message, name, level) \
{ \
  if (feasst::VERBOSE_LEVEL >= level) { \
    std::stringstream ss; \
    ss << " [" << feasst::feasst_dir_trim_(__FILE__) << ":" << __LINE__ \
              << "] " << message; \
    feasst_macro_output(name, ss.str()); \
  } \
}

void feasst_macro_output(const std::string& name, std::string message);

/// Define and name the various levels of verbosity.
# define WARN(message) FEASST_MACRO_OUTPUT(message, "Warn ", 2)
# define INFO(message) FEASST_MACRO_OUTPUT(message, "Info ", 3)
# define DEBUG(message) FEASST_MACRO_OUTPUT(message, "Debug", 4)
# define TRACE(message) FEASST_MACRO_OUTPUT(message, "Trace", 5)
}  // namespace feasst

#endif  // FEASST_UTILS_DEBUG_H_
