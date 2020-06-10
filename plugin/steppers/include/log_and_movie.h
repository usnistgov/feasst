
#ifndef FEASST_STEPPERS_LOG_AND_MOVIE_H_
#define FEASST_STEPPERS_LOG_AND_MOVIE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/analyze_factory.h"

namespace feasst {

/// Combine Log and Movie for convenience.
class LogAndMovie : public AnalyzeFactory {
 public:
  /**
    args:
    - file_name: .txt is appended to Log and .xyz is appended to Movie.
   */
  explicit LogAndMovie(const argtype& args = argtype());
};

inline std::shared_ptr<LogAndMovie> MakeLogAndMovie(
    const argtype &args = argtype()) {
  return std::make_shared<LogAndMovie>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_LOG_AND_MOVIE_H_
