
#ifndef FEASST_MATH_RANDOM_MT19937_H_
#define FEASST_MATH_RANDOM_MT19937_H_

#include "math/include/random.h"

namespace feasst {

/**
  Mersenne Twister 19937 generator
  http://www.cplusplus.com/reference/random/mt19937/
 */
class RandomMT19937 : public Random {
 public:
  RandomMT19937() { reseed_(); }

  // copy constructor should re-seed
  RandomMT19937(const RandomMT19937& random) { reseed_(); }

  /// Return a random real number with a uniform probability distribution
  /// between 0 and 1.
  double uniform() override { return dis_double_(generator_); }

  /// Return a random integer with a uniform probability distribution
  /// betwee min and max.
  int uniform(const int min, const int max) override;

  std::shared_ptr<Random> create(std::istream& istr) const override {
    return std::make_shared<RandomMT19937>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit RandomMT19937(std::istream& istr);
  virtual ~RandomMT19937() {}

 private:
  const std::string class_name_ = "RandomMT19937";
  std::uniform_real_distribution<double> dis_double_;
  std::mt19937 generator_;

  void reseed_() override;
};

inline std::shared_ptr<RandomMT19937> MakeRandomMT19937() {
  return std::make_shared<RandomMT19937>();
}

}  // namespace feasst

#endif  // FEASST_MATH_RANDOM_MT19937_H_
