
#ifndef FEASST_SERVER_SERVER_H_
#define FEASST_SERVER_SERVER_H_

#include <vector>
#include <memory>

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
 */
class Server {
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

  void serialize(std::ostream& ostr) const;
  explicit Server(std::istream& istr);
  ~Server();

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

}  // namespace feasst

#endif  // FEASST_SERVER_SERVER_H_
