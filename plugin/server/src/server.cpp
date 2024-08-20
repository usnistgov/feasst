#include <fstream>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <string.h>
#include <iostream>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "server/include/server.h"

//using namespace std;

namespace feasst {

Server::Server(argtype * args) {
  port_ = integer("port", args, 54321);
  buffer_size_ = integer("buffer_size", args, 1000);
  ASSERT(buffer_size_ % 2 == 0, "buffer_size must be even.");
}
Server::Server(argtype args) : Server(&args) {
  feasst_check_all_used(args);
}

Server::~Server() {
  if (bound_) {
    disconnect_();
  }
}

Server::Server(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 6579 && version <= 6579, "mismatch version: " << version);
  feasst_deserialize(&port_, istr);
  feasst_deserialize(&buffer_size_, istr);
}

void Server::serialize(std::ostream& ostr) const {
  feasst_serialize_version(6579, ostr);
  feasst_serialize(port_, ostr);
  feasst_serialize(buffer_size_, ostr);
}

// implementation from https://stackoverflow.com/questions/20732980/how-to-use-socket-with-a-python-client-and-a-c-server
void Server::bind_listen_accept() {
  if (bound_) {
    disconnect_();
  }
  std::cout << "# initializing server on localhost:" << port() << std::endl;
  buffer_ = new char[buffer_size() + 1];
  server_socket_= socket(AF_INET, SOCK_STREAM, 0);
  sockaddr_in server_addr;
  server_addr.sin_family = AF_INET;
  server_addr.sin_port = htons(port());
  server_addr.sin_addr.s_addr = INADDR_ANY;
  bind(server_socket_, (struct sockaddr*)&server_addr, sizeof(struct sockaddr));
  listen(server_socket_, 1);
  sockaddr_in clientAddr;
  socklen_t sin_size = sizeof(struct sockaddr_in);
  client_socket_ = accept(server_socket_, (struct sockaddr*)&clientAddr, &sin_size);
  bound_ = true;
//  // obtain ip address
//  {
//    struct sockaddr_in* pV4Addr = (struct sockaddr_in*)&serverAddr;
//    struct in_addr ipAddr = pV4Addr->sin_addr;
//    char ipstr[INET_ADDRSTRLEN];
//    inet_ntop( AF_INET, &ipAddr, ipstr, INET_ADDRSTRLEN );
//    INFO("ip " << ipstr);
//  }
}

int Server::receive() {
  bzero(buffer_, buffer_size_);
  const int size = read(client_socket_, buffer_, buffer_size_/2);
  DEBUG("received: " << buffer_);
  DEBUG("receive size: " << size);
  ASSERT(size >= 0, "error");
  return size;
}

void Server::send(const char* message) {
  strcpy(buffer_, message);
  const int size = write(client_socket_, buffer_, strlen(buffer_));
  DEBUG("sent: " << buffer_);
  DEBUG("send size: " << size);
}

void Server::send(const std::string message) {
  strcpy(buffer_, message.c_str());
  const size_t len = strlen(buffer_);
  ASSERT(len <= message.length(), "message is larger than buffer");
  const int size = write(client_socket_, buffer_, strlen(buffer_));
  DEBUG("sent: " << buffer_);
  DEBUG("send size: " << size);
}

//void Server::send_buffer() {
//  const int size = write(client_socket_, buffer_, strlen(buffer_));
//  DEBUG("sent: " << buffer_);
//  DEBUG("send size: " << size);
//}

void Server::disconnect_() {
  std::cout << "# closing server on localhost:" << port() << std::endl;
  close(server_socket_);
  delete[] buffer_;
}

}  // namespace feasst
