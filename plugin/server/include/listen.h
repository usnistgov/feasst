
#ifndef FEASST_SERVER_LISTEN_H_
#define FEASST_SERVER_LISTEN_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/action.h"

namespace feasst {

class Server;

/**
  Client server interface.
  Listen to receive a message from client.
  Each message should be processed as a single line with the same syntax as the
  text interface.
  This Action is not recommended for use in the text interface.
  Instead, use: echo "Server port 54321 buffer_size 1000" | /path/to/feasst/build/bin/fst
 */
class Listen : public Action {
 public:
  //@{
  /** @name Arguments
    - Server arguments.
   */
  explicit Listen(argtype args = argtype());
  explicit Listen(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<Listen>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<Listen>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Listen(std::istream& istr);
  virtual ~Listen() {}

  //@}
 private:
  std::shared_ptr<Server> server_;
};

inline std::shared_ptr<Listen> MakeListen(argtype args = argtype()) {
  return std::make_shared<Listen>(args);
}

}  // namespace feasst

#endif  // FEASST_SERVER_LISTEN_H_
