
#ifndef FEASST_MATH_RANDOM_MT19937_H_
#define FEASST_MATH_RANDOM_MT19937_H_

#include "math/include/random.h"

namespace feasst {

/**
  Mersenne Twister 19937 generator.
  http://www.cplusplus.com/reference/random/mt19937/
 */
class RandomMT19937 : public Random {
 public:
  RandomMT19937(const argtype& args = argtype());

  // serialize
  std::shared_ptr<Random> create(std::istream& istr) const override {
    return std::make_shared<RandomMT19937>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit RandomMT19937(std::istream& istr);
  virtual ~RandomMT19937() {}

 private:
  std::uniform_real_distribution<double> dis_double_;
  std::mt19937 generator_;

  void reseed_() override;

  double gen_uniform_() override { return dis_double_(generator_); }

  int gen_uniform_(const int min, const int max) override;
};

inline std::shared_ptr<RandomMT19937> MakeRandomMT19937(
    const argtype& args = argtype()) {
  return std::make_shared<RandomMT19937>(args);
}

}  // namespace feasst

#endif  // FEASST_MATH_RANDOM_MT19937_H_
