
#ifndef FEASST_SERVER_SERVER_H_
#define FEASST_SERVER_SERVER_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/action.h"

namespace feasst {

/**
 */
class Server : public Action {
 public:
  //@{
  /** @name Arguments
    - port: server port to listen for client (default: 54321).
    - buffer_size: maximum number of characters to fill message buffer.
      Must be even (default: 1000).
   */
  explicit Server(argtype args = argtype());
  explicit Server(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  int port() const { return port_; }
  int buffer_size() const { return buffer_size_; }
  const char* buffer() const { return buffer_; }
  void bind_listen_accept();
  bool bound() const { return bound_; }
  int receive();
  void send(const std::string message);
  void send(const char* message);
  //char * get_buffer() { return buffer_; }
  //void send_buffer();

  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<Server>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<Server>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Server(std::istream& istr);
  virtual ~Server();

  //@}
 private:
  int port_;
  int buffer_size_;
  int server_socket_;
  int client_socket_;
  char* buffer_;

  bool bound_ = false;
  void disconnect_();
};

inline std::shared_ptr<Server> MakeServer(argtype args = argtype()) {
  return std::make_shared<Server>(args);
}

}  // namespace feasst

#endif  // FEASST_SERVER_SERVER_H_
