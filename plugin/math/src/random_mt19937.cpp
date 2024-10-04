#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/random_mt19937.h"

namespace feasst {

FEASST_MAPPER(RandomMT19937, argtype());

RandomMT19937::RandomMT19937(argtype * args) : Random(args) {
  class_name_ = "RandomMT19937";
  std_normal_ = std::normal_distribution<double>(0., 1.);
  parse_seed_(args);
}
RandomMT19937::RandomMT19937(argtype args) : RandomMT19937(&args) {
  feasst_check_all_used(args);
}

void RandomMT19937::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_random_(ostr);
  feasst_serialize_version(101, ostr);
  ostr << MAX_PRECISION;
  ostr << generator_ << " ";
  feasst_serialize_endcap("RandomMT19937", ostr);
}

RandomMT19937::RandomMT19937(std::istream& istr)
  : Random(istr) {
  ASSERT(class_name_ == "RandomMT19937", "name: " << class_name_);
  const int verison = feasst_deserialize_version(istr);
  ASSERT(verison == 101, "version");
  istr >> generator_;
  feasst_deserialize_endcap("RandomMT19937", istr);
}

int RandomMT19937::gen_uniform_(const int min, const int max) {
  auto dis_int = std::uniform_int_distribution<int>(min, max);
  return dis_int(generator_);
}

void RandomMT19937::reseed_(const int seed) {
//  const int seed = rand();
  TRACE("seed " << seed << " address " << this);
  generator_ = std::mt19937(seed);
  dis_double_ = std::uniform_real_distribution<double>(0.0, 1.0);
}

}  // namespace feasst
