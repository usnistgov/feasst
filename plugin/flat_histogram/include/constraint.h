
#ifndef FEASST_FLAT_HISTOGRAM_CONSTRAINT_H_
#define FEASST_FLAT_HISTOGRAM_CONSTRAINT_H_

#include <memory>
#include <string>
#include "system/include/system.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/criteria.h"

namespace feasst {

// HWH consider moving constraints to MonteCarlo
/**
  Impose constraints on the System and Criteria.
 */
class Constraint {
 public:
  Constraint() {}
  virtual bool is_allowed(const System* system,
    const Criteria* criteria,
    const Acceptance& acceptance) const = 0;

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Constraint> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<Constraint> >& deserialize_map();
  std::shared_ptr<Constraint> deserialize(std::istream& istr);
  virtual ~Constraint() {}

 protected:
  std::string class_name_;
  void serialize_constraint_(std::ostream& ostr) const;
  explicit Constraint(std::istream& istr);
};

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_CONSTRAINT_H_
