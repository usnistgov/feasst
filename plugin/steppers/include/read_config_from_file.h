
#ifndef FEASST_STEPPERS_READ_CONFIG_FROM_FILE_H_
#define FEASST_STEPPERS_READ_CONFIG_FROM_FILE_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  This class is used for post processing a configuration file to analyze with a
  Stepper.
  Read through a number of configurations in a file.
  For each update, set the configuration to the next.
 */
class ReadConfigFromFile : public ModifyUpdateOnly {
 public:
  explicit ReadConfigFromFile(argtype args = argtype());
  explicit ReadConfigFromFile(argtype * args);
  void update(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override { criteria->update(); }

  // serialize
  std::string class_name() const override {
    return std::string("ReadConfigFromFile"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<ReadConfigFromFile>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<ReadConfigFromFile>(args); }
  ReadConfigFromFile(std::istream& istr);
};

inline std::shared_ptr<ReadConfigFromFile> MakeReadConfigFromFile(
    argtype args = argtype()) {
  return std::make_shared<ReadConfigFromFile>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_READ_CONFIG_FROM_FILE_H_
