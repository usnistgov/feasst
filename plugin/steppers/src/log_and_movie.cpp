#include "utils/include/serialize.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/log_and_movie.h"

namespace feasst {

LogAndMovie::LogAndMovie(argtype args) : AnalyzeFactory() {
  argtype log_args = args;
  feasst::append("file_name", &log_args, ".txt");
  add(MakeLog(log_args));
  feasst::append("file_name", &args, ".xyz");
  add(MakeMovie(args));
  WARN("LogAndMovie is deprecated. Please use Log and Movie separately.");
}

}  // namespace feasst
