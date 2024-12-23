
#ifndef FEASST_UTILS_INCLUDE_CHECKPOINT_H_
#define FEASST_UTILS_INCLUDE_CHECKPOINT_H_

#include <fstream>
#include <string>
#include <memory>
#include <map>
#include "utils/include/timer.h"
#include "utils/include/file.h"
#include "utils/include/debug.h"
#include "utils/include/io.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Save the state of a class in memory by writing to disk, such that the
  checkpoint file can be later read to restart the simulation.
  Note that for OMP or parallel simulations, the number of hours is multiplied
  by the number of threads.

  To restart a simulation from a checkpoint file, checkpoint.fst, use the
  following BASH command:

  echo "Restart checkpoint.fst" | $HOME/feasst(-version)/build/bin/fst

  Note that more commands may be added afted the line "Restart checkpoint.fst"
 */
class Checkpoint {
 public:
  //@{
  /** @name Arguments
    - num_hours: Number of hours between printing of checkpoint file
      (default: 1).
    - num_hours_terminate: Terminate after this many hours. If -1, do not
      terminate (default: -1).
      Termination may be detected in Bash shell using "$? != 0".
    - checkpoint_file: The default is one space (e.g., " ").
      Do not checkpoint if checkpoint_file is empty or is one space.
    - writes_per_backup: Create a unique checkpoint file name every this many
      times that a checkpoint is written (default: -1).
      If -1, only backup the previous file by appending its name with ".bak".
      Otherwise, if > 0, append each backup with an integer count beginning
      with 0.
   */
  explicit Checkpoint(argtype args = argtype());

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return number of hours between writing file.
  double num_hours() const { return num_hours_; }

  /// Write the checkpoint to file. If the file exists, create backup.
  template <typename T>
  void write(const T& obj, const std::string append_backup = ".bak") const {
    if (checkpoint_file_.empty() || checkpoint_file_ == " ") return;
    file_backup(checkpoint_file_, append_backup);
    std::ofstream file(checkpoint_file_.c_str(),
      std::ofstream::out | std::ofstream::trunc);
    std::stringstream ss;
    obj.serialize(ss);
    file << ss.str();
    file.close();
  }

  /// Write object to checkpoint_file if num_hours has passed since previous.
  template <typename T>
  void check(const T& obj) {
    const double hours = cpu_hours();
    const bool is_write = hours > previous_hours_ + num_hours_;
    const bool is_terminate = num_hours_terminate_ > 0 &&
                              hours > first_hours_ + num_hours_terminate_;
    if (is_write) previous_hours_ = hours;
    if (is_write || is_terminate) {
      std::string append_backup = ".bak";
      if (writes_per_backup_ > 0) {
        ++previous_backup_;
        file_backup(checkpoint_file_, feasst::str(previous_backup_));
      }
      write(obj, append_backup);
    }
    if (is_terminate) {
      FATAL("Terminating because Checkpoint has reached the user input " <<
        "num_hours_terminate: " << num_hours_terminate_ << ". Detect this " <<
        "termination in Bash shell using \"$? != 0\"");
    }
  }

  /// Initialize object by reading from file.
  template <typename T>
  void read(T * obj) {
    std::ifstream file(checkpoint_file_.c_str());
    ASSERT(file.good(), "cannot find " << checkpoint_file_);
    std::string line;
    std::getline(file, line);
    std::stringstream ss(line);
    *obj = T(ss);
  }
  template <typename T>
  void read_unique(std::unique_ptr<T>& obj) {
    std::ifstream file(checkpoint_file_.c_str());
    ASSERT(file.good(), "cannot find " << checkpoint_file_);
    std::string line;
    std::getline(file, line);
    std::stringstream ss(line);
    obj = std::make_unique<T>(ss);
  }

  /// Serialize object.
  void serialize(std::ostream& ostr) const;

  /// Deserialize object.
  explicit Checkpoint(std::istream& istr);

  //@}
 private:
  std::string checkpoint_file_;
  double num_hours_ = 0;
  double num_hours_terminate_ = 0;
  int writes_per_backup_;
  int previous_backup_ = -1;

  // temporary, not to be checkpointed
  double first_hours_ = -1.;
  double previous_hours_ = 0.;
};

inline std::shared_ptr<Checkpoint> MakeCheckpoint(argtype args = argtype()) {
  return std::make_shared<Checkpoint>(args);
}

}  // namespace feasst

#endif  // FEASST_UTILS_INCLUDE_CHECKPOINT_H_
