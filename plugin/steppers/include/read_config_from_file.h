
#ifndef FEASST_STEPPERS_READ_CONFIG_FROM_FILE_H_
#define FEASST_STEPPERS_READ_CONFIG_FROM_FILE_H_

#include <string>
#include <fstream>
#include "configuration/include/file_xyz.h"
#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  This class is used for post processing a configuration file to analyze with a
  Stepper.
  Read through a number of configurations in a file.
  For each update, set the configuration to the next.
  Once the end of file is reached, the Criteria is set to complete.
  Thus, use with "Run until_criteria_complete true"
 */
class ReadConfigFromFile : public ModifyUpdateOnly {
 public:
  //@{
  /** @name Arguments
    - input_file: name of FileXYZ to input Configuration.
    - Stepper arguments.
   */
  explicit ReadConfigFromFile(argtype args = argtype());
  explicit ReadConfigFromFile(argtype * args);
  //@}
  /** @name Public Functions
   */
  //@{
  void initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;
  void update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) override;

  // serialize
  std::string class_name() const override {
    return std::string("ReadConfigFromFile"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<ReadConfigFromFile>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<ReadConfigFromFile>(args); }
  ReadConfigFromFile(std::istream& istr);

  //@}
 private:
  std::string input_file_;
  bool set_complete_next_update_ = false;
  FileXYZ xyz_;

  // not serialized
  std::ifstream file_;

  void load_(Criteria * criteria, System * system);
};

inline std::shared_ptr<ReadConfigFromFile> MakeReadConfigFromFile(
    argtype args = argtype()) {
  return std::make_shared<ReadConfigFromFile>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_READ_CONFIG_FROM_FILE_H_
