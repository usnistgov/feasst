
#ifndef FEASST_MONTE_CARLO_CONSTRAINT_H_
#define FEASST_MONTE_CARLO_CONSTRAINT_H_

#include <memory>
#include <string>
#include "system/include/system.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/criteria.h"

namespace feasst {

/**
  Impose constraints on the System and Criteria.
 */
class Constraint {
 public:
  Constraint() {}
  virtual bool is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const = 0;

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Constraint> create(std::istream& istr) const;
  virtual std::shared_ptr<Constraint> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<Constraint> >& deserialize_map();
  std::shared_ptr<Constraint> deserialize(std::istream& istr);
  std::shared_ptr<Constraint> factory(const std::string name, argtype * args);
  virtual ~Constraint() {}

 protected:
  std::string class_name_;
  void serialize_constraint_(std::ostream& ostr) const;
  explicit Constraint(std::istream& istr);
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_CONSTRAINT_H_
