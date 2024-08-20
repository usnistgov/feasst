
#ifndef FEASST_MONTE_CARLO_ACTION_H_
#define FEASST_MONTE_CARLO_ACTION_H_

#include <memory>
#include <map>
#include <string>

namespace feasst {

class MonteCarlo;

typedef std::map<std::string, std::string> argtype;

/**
  An Action is for use with the MonteCarlo(arglist) constructor.
  Actions are performed immediately when they reach the top of the arglist.
  Examples include modifications to the simulations, or the initiation of a run.
 */
class Action {
 public:
  explicit Action(argtype args = argtype());
  explicit Action(argtype * args);
  virtual void run(MonteCarlo * mc);
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<Action> create(std::istream& istr) const;
  virtual std::shared_ptr<Action> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<Action> >& deserialize_map();
  std::shared_ptr<Action> deserialize(std::istream& istr);
  std::shared_ptr<Action> factory(const std::string name, argtype * args);
  virtual ~Action() {}

 protected:
  std::string class_name_ = "Action";
  void serialize_action_(std::ostream& ostr) const;
  explicit Action(std::istream& istr);
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_ACTION_H_
