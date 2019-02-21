#include <string>
#include <sstream>
#include "core/include/random.h"
#include "core/include/debug.h"

namespace feasst {

void seed_random_by_date() {
  const int t = time(NULL);
  srand ( t );
  INFO("time(seed): " << t);
}

void seed_random(const int seed) {
  srand ( seed );
  INFO("Initializing random number generator for reproduction with seed("
    << seed << ")");
}

int Random::uniform(const int min, const int max) {
  auto dis_int = std::uniform_int_distribution<int>(min, max);
  return dis_int(generator_);
}

// acknowledge: http://stackoverflow.com/questions/440133/
// how-do-i-create-a-random-alpha-numeric-string-in-c
std::string Random::alpha_numeric(const int length) {
  std::string str;
  const std::string alphanum = "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  for (int index = 0; index < length; ++index) {
    std::stringstream ss;
    ss << alphanum[uniform(0, alphanum.size() - 1)];
    str.append(ss.str());
  }
  return str;
}

double Random::uniform_real(const double min, const double max) {
  return (max - min)*uniform() + min;
}

}  // namespace feasst
