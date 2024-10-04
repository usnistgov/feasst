#include <fstream>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <string.h>
#include "utils/include/arguments_extra.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "monte_carlo/include/monte_carlo.h"
#include "server/include/server.h"
#include "server/include/listen.h"

//using namespace std;

namespace feasst {

Listen::Listen(argtype * args) {
  class_name_ = "Listen";
  server_ = std::make_unique<Server>(args);
}
Listen::Listen(argtype args) : Listen(&args) {
  feasst_check_all_used(args);
}
Listen::~Listen() {}

FEASST_MAPPER(Listen,);

Listen::Listen(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 3204 && version <= 3204, "mismatch version: " << version);
  feasst_deserialize(server_, istr);
}

void Listen::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3204, ostr);
  feasst_serialize(server_, ostr);
}

void Listen::run(MonteCarlo * mc) {
  server_->bind_listen_accept();
  bool finished = false;
  while (!finished) {
    const int size = server_->receive();
    if (strcmp(server_->buffer(), "EndListen") == 0) {
      DEBUG(server_->buffer());
      finished = true;
    } else {
      if (size > 0) {
        std::string line(server_->buffer());
        std::pair<std::string, argtype> marg = parse_line(line, NULL, NULL);
        DEBUG(marg.first);
        DEBUG(str(marg.second));
        arglist list(1);
        list[0] = marg;
        mc->parse_args(&list);
        ASSERT(list.size() == 0, "Unrecognized argument: " << list.begin()->first);
      }
      server_->send("Message received. Continue.");
    }
    ASSERT(size >= 0, "error");
  }
}

}  // namespace feasst
