#include "utils/include/serialize.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/log_and_movie.h"

namespace feasst {

LogAndMovie::LogAndMovie(const argtype& args) : AnalyzeFactory() {
  add(MakeLog(Arguments().append(".txt", "file_name", args)));
  add(MakeMovie(Arguments().append(".xyz", "file_name", args)));
}

}  // namespace feasst
