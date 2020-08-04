#include "math/include/random_mt19937.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"

namespace feasst {

int RandomMT19937::gen_uniform_(const int min, const int max) {
  auto dis_int = std::uniform_int_distribution<int>(min, max);
  return dis_int(generator_);
}

void RandomMT19937::reseed_() {
  const int seed = rand();
  TRACE("seed " << seed << " address " << this);
  generator_ = std::mt19937(seed);
  dis_double_ = std::uniform_real_distribution<double>(0.0, 1.0);
}

class MapRandomMT19937 {
 public:
  MapRandomMT19937() {
    RandomMT19937().deserialize_map()["RandomMT19937"] = MakeRandomMT19937();
  }
};

static MapRandomMT19937 mapper_ = MapRandomMT19937();

RandomMT19937::RandomMT19937(const argtype& args) : Random(args) {
  class_name_ = "RandomMT19937";
  parse_seed_(args);
}

void RandomMT19937::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_random_(ostr);
  feasst_serialize_version(101, ostr);
  ostr << MAX_PRECISION;
  ostr << generator_ << " ";
}

RandomMT19937::RandomMT19937(std::istream& istr)
  : Random(istr) {
  ASSERT(class_name_ == "RandomMT19937", "name: " << class_name_);
  const int verison = feasst_deserialize_version(istr);
  ASSERT(verison == 101, "version");
  istr >> generator_;
}

}  // namespace feasst
